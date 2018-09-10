// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights
// reserved. See files LICENSE and NOTICE for details.
//
// This file is part of CEED, a collection of benchmarks, miniapps, software
// libraries and APIs for efficient high-order finite element and spectral
// element discretizations for exascale applications. For more information and
// source code availability see http://github.com/ceed.
//
// The CEED research is supported by the Exascale Computing Project 17-SC-20-SC,
// a collaborative effort of two U.S. Department of Energy organizations (Office
// of Science and the National Nuclear Security Administration) responsible for
// the planning and preparation of a capable exascale ecosystem, including
// software, applications, hardware, advanced system engineering and early
// testbed platforms, in support of the nation's exascale computing imperative.

#include "c7_solver.hpp"

#ifdef MFEM_USE_MPI

using namespace std;

namespace mfem
{

namespace nth
{

void VisualizeField(socketstream &sock, const char *vishost, int visport,
                    ParGridFunction &gf, const char *title,
                    int x, int y, int w, int h, bool vec)
{
   ParMesh &pmesh = *gf.ParFESpace()->GetParMesh();
   MPI_Comm comm = pmesh.GetComm();

   int num_procs, myid;
   MPI_Comm_size(comm, &num_procs);
   MPI_Comm_rank(comm, &myid);

   bool newly_opened = false;
   int connection_failed;

   do
   {
      if (myid == 0)
      {
         if (!sock.is_open() || !sock)
         {
            sock.open(vishost, visport);
            sock.precision(8);
            newly_opened = true;
         }
         sock << "solution\n";
      }

      pmesh.PrintAsOne(sock);
      gf.SaveAsOne(sock);

      if (myid == 0 && newly_opened)
      {
         sock << "window_title '" << title << "'\n"
              << "window_geometry "
              << x << " " << y << " " << w << " " << h << "\n"
              << "keys maaAcl";
         if ( vec ) { sock << "vvv"; }
         sock << endl;
      }

      if (myid == 0)
      {
         connection_failed = !sock && !newly_opened;
      }
      MPI_Bcast(&connection_failed, 1, MPI_INT, 0, comm);
   }
   while (connection_failed);
}

C7Operator::C7Operator(int size, 
                       ParFiniteElementSpace &h1_fes,
                       ParFiniteElementSpace &l2_fes, 
                       Array<int> &essential_tdofs,
                       ParGridFunction &rho0, 
                       double cfl_, 
                       AWBSMasterOfPhysics *AWBSPhysics_,
					   ParGridFunction &x_gf_, 
                       ParGridFunction &T_gf_, 
                       double cgt, int cgiter)
   : TimeDependentOperator(size),
     H1FESpace(h1_fes), L2FESpace(l2_fes),
     H1compFESpace(h1_fes.GetParMesh(), h1_fes.FEColl(), 1),
     ess_tdofs(essential_tdofs),
     dim(h1_fes.GetMesh()->Dimension()),
     nzones(h1_fes.GetMesh()->GetNE()),
     l2dofs_cnt(l2_fes.GetFE(0)->GetDof()),
     h1dofs_cnt(h1_fes.GetFE(0)->GetDof()),
     cfl(cfl_), cg_rel_tol(cgt), cg_max_iter(cgiter), M1_closure(false),
     invM0nu(&l2_fes), M1nu(&h1_fes), M1nut(&h1_fes), B1(&h1_fes), 
     DI(&l2_fes, &h1_fes), VEscaled(&l2_fes, &h1_fes), 
	 VEscatter(&l2_fes, &h1_fes), 
     DA(&l2_fes, &h1_fes), VAEscaled(&l2_fes, &h1_fes), 
     VAEfull(&l2_fes, &h1_fes),
     integ_rule(IntRules.Get(h1_fes.GetMesh()->GetElementBaseGeometry(),
                             3*h1_fes.GetOrder(0) + l2_fes.GetOrder(0) - 1)),
     quad_data(dim, nzones, integ_rule.GetNPoints()),
     quad_data_is_current(false),
     timer(), AWBSPhysics(AWBSPhysics_), x_gf(x_gf_)
{
   GridFunctionCoefficient rho_coeff(&rho0);

   // Values of rho0DetJ0 and Jac0inv at all quadrature points.
   const int nqp = integ_rule.GetNPoints();
   Vector rho_vals(nqp);
   for (int i = 0; i < nzones; i++)
   {
      rho0.GetValues(i, integ_rule, rho_vals);
      ElementTransformation *T = h1_fes.GetElementTransformation(i);
      for (int q = 0; q < nqp; q++)
      {
         const IntegrationPoint &ip = integ_rule.IntPoint(q);
         T->SetIntPoint(&ip);

         DenseMatrixInverse Jinv(T->Jacobian());
         Jinv.GetInverseMatrix(quad_data.Jac0inv(i*nqp + q));

         const double rho0DetJ0 = T->Weight() * rho_vals(q);
         quad_data.rho0DetJ0w(i*nqp + q) = rho0DetJ0 *
                                           integ_rule.IntPoint(q).weight;
      }
   }

   // Initial local mesh size (assumes similar cells).
   double loc_area = 0.0, glob_area;
   int loc_z_cnt = nzones, glob_z_cnt;
   ParMesh *pm = H1FESpace.GetParMesh();
   for (int i = 0; i < nzones; i++) { loc_area += pm->GetElementVolume(i); }
   MPI_Allreduce(&loc_area, &glob_area, 1, MPI_DOUBLE, MPI_SUM, pm->GetComm());
   MPI_Allreduce(&loc_z_cnt, &glob_z_cnt, 1, MPI_INT, MPI_SUM, pm->GetComm());
   switch (pm->GetElementBaseGeometry(0))
   {
      case Geometry::SEGMENT:
         quad_data.h0 = glob_area / glob_z_cnt; break;
      case Geometry::SQUARE:
         quad_data.h0 = sqrt(glob_area / glob_z_cnt); break;
      case Geometry::TRIANGLE:
         quad_data.h0 = sqrt(2.0 * glob_area / glob_z_cnt); break;
      case Geometry::CUBE:
         quad_data.h0 = pow(glob_area / glob_z_cnt, 1.0/3.0); break;
      case Geometry::TETRAHEDRON:
         quad_data.h0 = pow(6.0 * glob_area / glob_z_cnt, 1.0/3.0); break;
      default: MFEM_ABORT("Unknown zone type!");
   }
   quad_data.h0 /= (double) H1FESpace.GetOrder(0);

   // Inverted Mf0 intergrators.
   invMass0NuIntegrator *invmf0nui = new invMass0NuIntegrator(quad_data);
   invmf0nui->SetIntRule(&integ_rule);
   invM0nu.AddDomainIntegrator(invmf0nui);

   // Standard assembly for the velocity mass matrix.
   Mass1NuIntegrator *mf1nui = new Mass1NuIntegrator(quad_data);
   mf1nui->SetIntRule(&integ_rule);
   M1nu.AddDomainIntegrator(mf1nui);
   M1nu.Assemble(); 

   Mass1NutIntegrator *mf1nuti = new Mass1NutIntegrator(quad_data);
   mf1nuti->SetIntRule(&integ_rule);
   M1nut.AddDomainIntegrator(mf1nuti);
   M1nut.Assemble();

   BfieldIntegrator *f1bfi = new BfieldIntegrator(quad_data);
   f1bfi->SetIntRule(&integ_rule);
   B1.AddDomainIntegrator(f1bfi);
   B1.Assemble();

   Divf1Integrator *f1di = new Divf1Integrator(quad_data);
   f1di->SetIntRule(&integ_rule);
   DA.AddDomainIntegrator(f1di);
   // Make a dummy assembly to figure out the sparsity.
   DA.Assemble(0);
   DA.Finalize(0);

   Divf0Integrator *tfi = new Divf0Integrator(quad_data);
   tfi->SetIntRule(&integ_rule);
   DI.AddDomainIntegrator(tfi);
   // Make a dummy assembly to figure out the sparsity.
   DI.Assemble(0);
   DI.Finalize(0);

   EfieldScatterIntegrator *_f0ei = new EfieldScatterIntegrator(quad_data);
   _f0ei->SetIntRule(&integ_rule);
   VEscatter.AddDomainIntegrator(_f0ei);
   // Make a dummy assembly to figure out the sparsity.
   VEscatter.Assemble(0);
   VEscatter.Finalize(0);

   EfieldScaledIntegrator *f0ei = new EfieldScaledIntegrator(quad_data);
   f0ei->SetIntRule(&integ_rule);
   VEscaled.AddDomainIntegrator(f0ei);
   // Make a dummy assembly to figure out the sparsity.
   VEscaled.Assemble(0);
   VEscaled.Finalize(0);

   AEfieldScaledIntegrator *f1aei = new AEfieldScaledIntegrator(quad_data);
   f1aei->SetIntRule(&integ_rule);
   VAEscaled.AddDomainIntegrator(f1aei);
   // Make a dummy assembly to figure out the sparsity.
   VAEscaled.Assemble(0);
   VAEscaled.Finalize(0);

   AEfieldIntegrator *f1Maei = new AEfieldIntegrator(quad_data);
   f1Maei->SetIntRule(&integ_rule);
   VAEfull.AddDomainIntegrator(f1Maei);
   // Make a dummy assembly to figure out the sparsity.
   VAEfull.Assemble(0);
   VAEfull.Finalize(0);
}

void C7Operator::Mult(const Vector &F, Vector &dFdv) const
{
   /////////////////////////////////////////////////////////////////////////////
   // We solve the first moment model of AWBS steady state electron transport.//
   /////////////////////////////////////////////////////////////////////////////
   ////// f0 equation //////////////////////////////////////////////////////////
   //                                                                         //
   // M0(nu)*df0dv - 1/v*M0(E*f1/f0)*df0dv = DT(I)*f1 + 2/v^2*VT(E)*f1        //
   //                                        + M0(nu)*dfMdv                   //
   //                                                                         //
   ////// f1 equation //////////////////////////////////////////////////////////
   //                                                                         //
   // M1(nu)*df1dv - 1/v*V(AE)*df0dv = -D(I)*f0 + 1/v^2*V((3A-I)E)*f0         //
   //                                  + 1/v*B(B)*f1 + 1/v*M1(nut)*f1         //
   //                                                                         //
   /////////////////////////////////////////////////////////////////////////////

   dFdv = 0.0;

   const double velocity = GetTime(); 

   UpdateQuadratureData(velocity, F);
   const double N_x_vTmax = AWBSPhysics->GetVelocityScale();
   const double velocity_real = velocity * N_x_vTmax;

   AWBSPhysics->dfMdv_pcf->SetVelocityReal(velocity_real);
   ParGridFunction dfMdv_source(&L2FESpace);
   dfMdv_source.ProjectCoefficient(*(AWBSPhysics->dfMdv_pcf));

   // The monolithic BlockVector stores the unknown fields as follows:
   // - isotropic F0 (energy density)
   // - anisotropic F1 (flux density)

   const int VsizeL2 = L2FESpace.GetVSize();
   const int VsizeH1 = H1FESpace.GetVSize();
   int size;

   ParGridFunction F0, F1, qH_gf, jC_gf;
   Vector* sptr = (Vector*) &F;
   size = 0;
   F0.MakeRef(&L2FESpace, *sptr, size);
   size += VsizeL2;
   F1.MakeRef(&H1FESpace, *sptr, size);
   size += VsizeH1;
   qH_gf.MakeRef(&H1FESpace, *sptr, size);
   size += VsizeH1;
   jC_gf.MakeRef(&H1FESpace, *sptr, size);

   ParGridFunction dF0, dF1, dqH_gf, djC_gf;
   size = 0;
   dF0.MakeRef(&L2FESpace, dFdv, 0);
   size += VsizeL2;
   dF1.MakeRef(&H1FESpace, dFdv, VsizeL2);
   size += VsizeH1;
   dqH_gf.MakeRef(&H1FESpace, dFdv, size);
   size += VsizeH1;
   djC_gf.MakeRef(&H1FESpace, dFdv, size);

   //Mf0nu.Update();
   //invM0nuE.Update();
   M1nu.Update(); 
   M1nut.Update();
   B1.Update();
   DI = 0.0;
   VEscaled = 0.0;
   DA = 0.0; 
   VAEscaled = 0.0;

   timer.sw_force.Start(); 
   //Mf0nu.Assemble();
   //invM0nuE.Assemble();
   M1nu.Assemble(); 
   M1nut.Assemble();
   B1.Assemble();
   DI.Assemble();
   VEscaled.Assemble(0);   
   DA.Assemble();  
   VAEscaled.Assemble(0); 
   timer.sw_force.Stop();

   // Solve for df0dv. 
   ////// f0 equation //////////////////////////////////////////////////////////
   //                                                                         //
   // M0(nu)*df0dv - 1/v*M0(E*f1/f0)*df0dv = DT(I)*f1 + 2/v^2*VT(E)*f1        //
   //                                        + M0(nu)*dfMdv                   //
   //                                                                         //
   /////////////////////////////////////////////////////////////////////////////
   Vector F0_rhs(VsizeL2);
   timer.sw_force.Start();
   DI.MultTranspose(F1, F0_rhs);
   VEscaled.AddMultTranspose(F1, F0_rhs, 
                             2.0 / velocity_real / velocity_real);
   //Mf0nu.AddMult(dfMdv_source, F0_rhs); 
   // Compute dF0.
   //invM0nuE.Mult(F0_rhs, dF0);

   // Solve for df1dv.
   ////// f1 equation //////////////////////////////////////////////////////////
   //                                                                         //
   // M1(nu)*df1dv - 1/v*V(AE)*df0dv = -D(I)*f0 + 1/v^2*V((3A-I)E)*f0         //
   //                                  + 1/v*B(B)*f1 + 1/v*M1(nut)*f1         //
   //                                                                         //
   /////////////////////////////////////////////////////////////////////////////
   Vector F1_rhs(VsizeH1), B, X;
   timer.sw_force.Start();
   DA.Mult(F0, F1_rhs);
   F1_rhs.Neg(); 
   B1.AddMult(F1, F1_rhs, 1.0 / velocity_real);
   M1nut.AddMult(F1, F1_rhs, 1.0 / velocity_real);
   VAEscaled.AddMult(dF0, F1_rhs, 1.0 / velocity_real);
   // Compute dF1.
   timer.sw_force.Stop();
   timer.dof_tstep += H1FESpace.GlobalTrueVSize();
   HypreParMatrix A;
   dF1 = 0.0;
   M1nu.FormLinearSystem(ess_tdofs, dF1, F1_rhs, A, X, B);
   CGSolver cg(H1FESpace.GetParMesh()->GetComm());
   cg.SetOperator(A);
   cg.SetRelTol(cg_rel_tol); cg.SetAbsTol(0.0);
   cg.SetMaxIter(cg_max_iter);
   cg.SetPrintLevel(0);
   timer.sw_cgH1.Start();
   cg.Mult(B, X);
   timer.sw_cgH1.Stop();
   timer.H1dof_iter += cg.GetNumIterations() * H1compFESpace.GlobalTrueVSize();
   M1nu.RecoverFEMSolution(X, F1_rhs, dF1);

   quad_data_is_current = false;

   // Integrate heat flux.
   dqH_gf = F1;
   dqH_gf *= pow(velocity_real, 5.0);
   // Since the integration goes from vmax -> vmin, revert sign.
   dqH_gf *= -1.0;
   // Integrate current.
   djC_gf = F1;
   djC_gf *= pow(velocity_real, 3.0);
   // Since the integration goes from vmax -> vmin, revert sign.
   djC_gf *= -1.0;

   // c7_oper uses a more general formulation of velocity space with scaled
   // velocity magnitude from v_normalized in (0, 1) to v_real in (0, NxvTmax),
   // where the maximum velocity is an N times multiple of thermal velocity 
   // given for the maximum value of the temperature profile.
   // Consequently, dFdv needs to be NxvTmax multiplied in order to provide 
   // a "real" integration F = int(dFdv, dv).
   dFdv *= N_x_vTmax;
}

void C7Operator::ImplicitSolve(const double dv, const Vector &F, Vector &dFdv)
{
   /////////////////////////////////////////////////////////////////////////////
   // We solve the first moment model of AWBS steady state electron transport.//
   /////////////////////////////////////////////////////////////////////////////
   ////// f0 equation //////////////////////////////////////////////////////////
   //                                                                         //
   // M0(nu)*df0dv - 1/v*VT(E)*df1dv = DT(I)*f1 + 2/v^2*VT(E)*f1              //
   //                                  + M0(nu)*dfMdv                         //
   //                                                                         //
   ////// f1 equation //////////////////////////////////////////////////////////
   //                                                                         //
   // M1(nu)*df1dv - 1/v*V(AE)*df0dv = -D(I)*f0 + 1/v^2*V((3A-I)E)*f0         //
   //                                  + 1/v*B(B)*f1 + 1/v*M1(nut)*f1         //
   //                                                                         //
   /////////////////////////////////////////////////////////////////////////////

   dFdv = 0.0;

   const double velocity = GetTime(); 

   UpdateQuadratureData(velocity, F);

   const double N_x_vTmax = AWBSPhysics->GetVelocityScale();
   // Get the real velocity and velocity step.
   const double velocity_real = velocity * N_x_vTmax;
   const double dv_real = dv * N_x_vTmax; 
 
   ParGridFunction dfMdv_source(&L2FESpace), fM_source(&L2FESpace);
   AWBSPhysics->dfMdv_pcf->SetVelocityReal(velocity_real);
   dfMdv_source.ProjectCoefficient(*(AWBSPhysics->dfMdv_pcf));
   AWBSPhysics->fM_pcf->SetVelocityReal(velocity_real);
   fM_source.ProjectCoefficient(*(AWBSPhysics->fM_pcf)); 

   // The monolithic BlockVector stores the unknown fields as follows:
   // - isotropic F0 (energy density)
   // - anisotropic F1 (flux density)

   const int VsizeL2 = L2FESpace.GetVSize();
   const int VsizeH1 = H1FESpace.GetVSize();
   int size;

   ParGridFunction F0, F1, qH_gf, jC_gf, a0_gf, b0_gf, b1_gf;
   Vector* sptr = (Vector*) &F;
   size = 0;
   F0.MakeRef(&L2FESpace, *sptr, size);
   size += VsizeL2;
   F1.MakeRef(&H1FESpace, *sptr, size);
   size += VsizeH1;
   qH_gf.MakeRef(&H1FESpace, *sptr, size);
   size += VsizeH1;
   jC_gf.MakeRef(&H1FESpace, *sptr, size);
   size += VsizeH1;
   a0_gf.MakeRef(&L2FESpace, *sptr, size);
   size += VsizeL2;
   b0_gf.MakeRef(&H1FESpace, *sptr, size);
   size += VsizeH1;
   b1_gf.MakeRef(&H1FESpace, *sptr, size);

   ParGridFunction dF0, dF1, dqH_gf, djC_gf, da0_gf, db0_gf, db1_gf;
   size = 0;
   dF0.MakeRef(&L2FESpace, dFdv, size);
   size += VsizeL2;
   dF1.MakeRef(&H1FESpace, dFdv, size);
   size += VsizeH1;
   dqH_gf.MakeRef(&H1FESpace, dFdv, size);
   size += VsizeH1;
   djC_gf.MakeRef(&H1FESpace, dFdv, size);
   size += VsizeH1;
   da0_gf.MakeRef(&L2FESpace, dFdv, size);
   size += VsizeL2;
   db0_gf.MakeRef(&H1FESpace, dFdv, size);
   size += VsizeH1;
   db1_gf.MakeRef(&H1FESpace, dFdv, size);

   invM0nu.Update();
   M1nu.Update(); 
   M1nut.Update();
   B1.Update();
   DI = 0.0;
   VEscatter = 0.0;
   VEscaled = 0.0;
   DA = 0.0; 
   VAEscaled = 0.0;
   VAEfull = 0.0;

   timer.sw_force.Start(); 
   invM0nu.Assemble();
   invM0nu.Finalize();
   M1nu.Assemble(); 
   M1nu.Finalize();
   M1nut.Assemble();
   M1nut.Finalize();
   B1.Assemble();
   B1.Finalize();
   DI.Assemble();
   DI.Finalize();
   VEscatter.Assemble(0);   
   VEscatter.Finalize();
   VEscaled.Assemble(0);   
   VEscaled.Finalize();
   DA.Assemble();  
   DA.Finalize();
   VAEscaled.Assemble(0); 
   VAEscaled.Finalize();
   VAEfull.Assemble(0); 
   VAEfull.Finalize();
   timer.sw_force.Stop();

   // Force no effect of Bfield.
   B1.SpMat()   *= 0.0; 

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////// Embedded iteration scheme //////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
   // The system matrix is constructed from M1nu mass matrix.
   ParBilinearForm &A_f1 = M1nu; // dF1 system matrix. Initialized as M1nu.

   // Transposed matrices.
   // Full E field matrix.
   SparseMatrix *_tVE = Transpose(VEscatter.SpMat());
   // Diffusion plus directional effect of Efield matrices.
   SparseMatrix *_tDIVE = Transpose(DI.SpMat()); 
   _tDIVE->Add(-2.0 / velocity_real / velocity_real, *_tVE);
   // Scaled E field matrix.
   SparseMatrix *_tVEsc = Transpose(VEscaled.SpMat());
   // Proceed with matrix inversions.
   SparseMatrix *_invM0_tDIVE = mfem::Mult(invM0nu.SpMat(), *_tDIVE);
   SparseMatrix *_invM0_tVEsc = mfem::Mult(invM0nu.SpMat(), *_tVEsc);
   SparseMatrix *_DA_invM0_tDIVE = mfem::Mult(DA.SpMat(), *_invM0_tDIVE);
   SparseMatrix *_DA_invM0_tVEsc = mfem::Mult(DA.SpMat(), *_invM0_tVEsc);

   // Complete the system matrix A_f1 initialized by the operator M1nu.
   A_f1.SpMat().Add(-1.0 * dv_real / velocity_real, M1nut.SpMat());
   A_f1.SpMat().Add(-1.0 * dv_real / velocity_real, B1.SpMat());
   A_f1.SpMat().Add(dv_real * dv_real, *_DA_invM0_tDIVE);

   // Full and partial RHS vectors.
   Vector b1(VsizeH1), b0(VsizeL2), b1_n(VsizeH1), b1_k0(VsizeH1), 
          b1_k1(VsizeH1), b0_n(VsizeL2), b0_k(VsizeL2), b0_kplus1(VsizeL2);
   // Fill b1_n vector.
   DA.Mult(fM_source, b1_n);
   DA.AddMult(F0, b1_n);
   b1_n.Neg();
   VAEfull.AddMult(dfMdv_source, b1_n, 1.0 / velocity_real);
   M1nut.AddMult(F1, b1_n, 1.0 / velocity_real);
   B1.AddMult(F1, b1_n, 1.0 / velocity_real);
   _DA_invM0_tDIVE->AddMult(F1, b1_n, -1.0 * dv_real);
   // Fill b0_n vector.
   _invM0_tDIVE->Mult(F1, b0_n);

   // Unknown at level k = 0 equal zero.
   //Vector dF1(VsizeH1), dF0(VsizeL2);
   dF1 = 0.0;
   dF0 = 0.0;

   double delta_dF_norm = 1.0; 
   double dF1_norm = 0.0;
   double converg_lim = 0.001;
   int kiter_max = 25;
   int k;
   int PCGNumIter_dF1;
   double glob_dF1_norm, glob_dF0_norm;
   double glob_b1_n_norm, glob_b1_k0_norm, glob_b1_k1_norm;
   for (k = 1; k < kiter_max && delta_dF_norm > converg_lim; k++)
   {
      // Fill b1_k0 and b1_k1 vectors.
      b1_k0 = 0.0;
      VAEscaled.AddMult(dF0, b1_k0, 1.0 / velocity_real);
	  b1_k1 = 0.0;
      _DA_invM0_tVEsc->AddMult(dF1, b1_k1, -1.0 * dv_real / velocity_real);
      // Fill b0_k vector. To be filled before dF1^k+1 has been solved.
      b0_k = 0.0;
      _invM0_tVEsc->AddMult(dF1, b0_k, 1.0 / velocity_real); // correct dF1^k

      // Fill full F1 RHS vectors.
      b1 = b1_n;
      b1 += b1_k0;
	  b1 += b1_k1;
      // Run the HYPRE solver.
      timer.sw_force.Stop();
      timer.dof_tstep += H1FESpace.GlobalTrueVSize();
      Vector _B, _X;
      HypreParMatrix _A; 
      dF1 = 0.0;    
      A_f1.FormLinearSystem(ess_tdofs, dF1, b1, _A, _X, _B);
      bool verbose = false;
      HypreBoomerAMG amg_dF1(_A);
      HyprePCG pcg_dF1(_A);
      pcg_dF1.SetTol(cg_rel_tol);
      pcg_dF1.SetMaxIter(cg_max_iter);
      pcg_dF1.SetPrintLevel(verbose);
      pcg_dF1.SetPreconditioner(amg_dF1);
      amg_dF1.SetPrintLevel(verbose);
      // Solve for dF1.
      pcg_dF1.Mult(_B, _X);
      timer.sw_cgH1.Stop();
      // Number of PCG iterations.
      pcg_dF1.GetNumIterations(PCGNumIter_dF1);
      timer.H1dof_iter += PCGNumIter_dF1 * H1compFESpace.GlobalTrueVSize();    
      A_f1.RecoverFEMSolution(_X, b1, dF1);

      // Fill b0_k vector. Trying to be fill after dF1^k+1 has been solved.
      //b0_k = 0.0;
      //_invM0_tVEsc->AddMult(dF1, b0_k, 1.0 / velocity_real); // correct dF1^k
      // Fill b0_k+1 vector. To be filled after dF1^k+1 has been solved.
      b0_kplus1 = 0.0;
      _invM0_tDIVE->AddMult(dF1, b0_kplus1, dv_real); // dF1^k+1

      // Fill full F0 RHS vector.
      b0 = b0_n;
      b0 += b0_k;
      b0 += b0_kplus1;

      // Solve for dF0.
      dF0 = b0;

      // L2 norms will be used for convergence.
      double loc_dF1_norm = dF1.Norml2();
      MPI_Allreduce(&loc_dF1_norm, &glob_dF1_norm, 1, MPI_DOUBLE, MPI_SUM, 
	                H1FESpace.GetParMesh()->GetComm());
      double loc_dF0_norm = dF0.Norml2();
      MPI_Allreduce(&loc_dF0_norm, &glob_dF0_norm, 1, MPI_DOUBLE, MPI_SUM, 
	                H1FESpace.GetParMesh()->GetComm());
	  // Relative difference in the norm evolution.
	  delta_dF_norm = abs((dF1_norm - glob_dF1_norm) / glob_dF1_norm);
	  dF1_norm = glob_dF1_norm;
      double loc_b1_n_norm = b1_n.Norml2(), loc_b1_k0_norm = b1_k0.Norml2(),
             loc_b1_k1_norm = b1_k1.Norml2();
      MPI_Allreduce(&loc_b1_n_norm, &glob_b1_n_norm, 1, MPI_DOUBLE, MPI_SUM, 
	                H1FESpace.GetParMesh()->GetComm());
      MPI_Allreduce(&loc_b1_k0_norm, &glob_b1_k0_norm, 1, MPI_DOUBLE, MPI_SUM, 
	                H1FESpace.GetParMesh()->GetComm());
      MPI_Allreduce(&loc_b1_k1_norm, &glob_b1_k1_norm, 1, MPI_DOUBLE, MPI_SUM, 
	                H1FESpace.GetParMesh()->GetComm());
      // Output on root processor.
      if (0)
	  {
      if (H1FESpace.GetParMesh()->GetMyRank() == 0)
      { 
         cout << "kth iteration: " <<  k << endl << flush;
         cout << "HyprePCG_dF1(BoomerAMG) GetNumIterations: " 
              <<  PCGNumIter_dF1 << endl << flush;
         cout << "dF1 L2 norm: " 
              <<  glob_dF1_norm //* dv_real * pow(velocity_real, 5.0)
              << endl << flush;
         cout << "dF0 L2 norm: " 
              <<  glob_dF0_norm //* dv_real * pow(velocity_real, 5.0)
              << endl << flush;
         cout << "delta_dF_norm: " <<  delta_dF_norm << endl << flush;
         cout << "b1_n_norm: " << scientific
              <<  glob_b1_n_norm 
              << endl << flush;
         cout << "b1_k0_norm: " << scientific
              <<  glob_b1_k0_norm 
              << endl << flush;
         cout << "b1_k1_norm: " << scientific
              <<  glob_b1_k1_norm 
              << endl << flush;
      }
	  }
   }
   // Output on root processor.
   if (H1FESpace.GetParMesh()->GetMyRank() == 0)
   { 
      cout << "kth iteration: " <<  k << endl << flush;
      cout << "HyprePCG_dF1(BoomerAMG) GetNumIterations: " 
           <<  PCGNumIter_dF1 << endl << flush;
      cout << "dF1 L2 norm: " << scientific
           <<  glob_dF1_norm //* dv_real * pow(velocity_real, 5.0)
           << endl << flush;
      cout << "dF0 L2 norm: " << scientific
           <<  glob_dF0_norm //* dv_real //* pow(velocity_real, 5.0)
           << endl << flush;
      cout << "delta_dF_norm: " <<  delta_dF_norm << endl << flush;
      if (delta_dF_norm > converg_lim)
      {
         cout << "Embedded E.dFdv itaration has not converged!" 
              << endl << flush;
         cout << "b1_n_norm: " << scientific
              <<  glob_b1_n_norm 
              << endl << flush;
         cout << "b1_k0_norm: " << scientific
              <<  glob_b1_k0_norm 
              << endl << flush;
         cout << "b1_k1_norm: " << scientific
              <<  glob_b1_k1_norm 
              << endl << flush;
      } 
   }
   
   // Clean the buffer.
   delete _tVE;
   delete _tDIVE; 
   delete _tVEsc;
   delete _invM0_tDIVE;
   delete _invM0_tVEsc;
   delete _DA_invM0_tDIVE;
   delete _DA_invM0_tVEsc;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////// Embedded iteration scheme //////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

   quad_data_is_current = false;

   // Integrate heat flux.
   dqH_gf = F1;
   dqH_gf *= pow(velocity_real, 5.0);
   // Integrate current.
   djC_gf = F1;
   djC_gf *= pow(velocity_real, 3.0);

   // Generalized Ohm's law related grid functions.
   AWBSPhysics->P1a0_pcf->SetdFM(&dfMdv_source);
   AWBSPhysics->P1a0_pcf->SetdF0(&dF0);
   AWBSPhysics->P1a0_pcf->SetVelocityReal(velocity_real); 
   Coefficient &P1a0_cf = *(AWBSPhysics->P1a0_pcf);
   da0_gf.ProjectCoefficient(P1a0_cf);

   ////////////////////////////////////////////////////////////////////////////
   ///// The weak solution to 1/nut*grad(fM + delta f0) ///////////////////////
   ////////////////////////////////////////////////////////////////////////////
   ParGridFunction gradAf0(&H1FESpace), gradAf0_rhs(&H1FESpace);
   // Weak derivative.
   DA.Mult(fM_source, gradAf0_rhs);
   DA.AddMult(F0, gradAf0_rhs);
   gradAf0_rhs.Neg();
   // Initialize unknown vector.
   gradAf0 = 0.0;
   Vector B_, X_;
   HypreParMatrix A_;
   M1nut.FormLinearSystem(ess_tdofs, gradAf0, gradAf0_rhs, A_, X_, B_);
   bool verbose_gradAf0 = false;
   HypreBoomerAMG amg_gradAf0(A_);
   HyprePCG pcg_gradAf0(A_);
   pcg_gradAf0.SetTol(cg_rel_tol);
   pcg_gradAf0.SetMaxIter(cg_max_iter);
   pcg_gradAf0.SetPrintLevel(verbose_gradAf0);
   pcg_gradAf0.SetPreconditioner(amg_gradAf0);
   amg_gradAf0.SetPrintLevel(verbose_gradAf0);
   pcg_gradAf0.Mult(B_, X_);
   timer.sw_cgH1.Stop();
   int PCGNumIter_gradAf0;
   pcg_gradAf0.GetNumIterations(PCGNumIter_gradAf0);
   timer.H1dof_iter += PCGNumIter_gradAf0 * H1compFESpace.GlobalTrueVSize();
   if (H1FESpace.GetParMesh()->GetMyRank() == 0)
   { 
      cout << "HyprePCG_gradAf0(BoomerAMG) GetNumIterations: " 
           <<  PCGNumIter_gradAf0 << endl << flush;
   }    
   M1nut.RecoverFEMSolution(X_, gradAf0_rhs, gradAf0);
   ////////////////////////////////////////////////////////////////////////////
   ///// The weak solution to 1/nut*grad(fM + delta f0) ///////////////////////
   ////////////////////////////////////////////////////////////////////////////

   // Evaluate the P1b0 coefficient.
   AWBSPhysics->P1b0_pcf->SetF1(&gradAf0);
   //fM_source = 0.0;
   // Put together fM and delta f0.
   //fM_source += F0;
   //AWBSPhysics->P1b0_pcf->SetF0(&fM_source);
   //AWBSPhysics->P1b0_pcf->SetF0(&F0);
   AWBSPhysics->P1b0_pcf->SetdF1(&dF1);
   AWBSPhysics->P1b0_pcf->SetVelocityReal(velocity_real);
   VectorCoefficient &P1b0_cf = *(AWBSPhysics->P1b0_pcf);
   db0_gf.ProjectCoefficient(P1b0_cf);

   AWBSPhysics->P1b1_pcf->SetF1(&dF1);
   AWBSPhysics->P1b1_pcf->SetVelocityReal(velocity_real);
   VectorCoefficient &P1b1_cf = *(AWBSPhysics->P1b1_pcf);
   db1_gf.ProjectCoefficient(P1b1_cf);

   // c7_oper uses a more general formulation of velocity space with scaled
   // velocity magnitude from v_normalized in (0, 1) to v_real in (0, NxvTmax),
   // where the maximum velocity is an N times multiple of thermal velocity 
   // given for the maximum value of the temperature profile.
   // Consequently, dFdv needs to be NxvTmax multiplied in order to provide 
   // a "real" integration F = int(dFdv, dv).
   dFdv *= N_x_vTmax;
}

double C7Operator::GetVelocityStepEstimate(const Vector &S) const
{
   const double velocity = GetTime();
   UpdateQuadratureData(velocity, S);

   double glob_dt_est;
   MPI_Allreduce(&quad_data.dt_est, &glob_dt_est, 1, MPI_DOUBLE, MPI_MIN,
                 H1FESpace.GetParMesh()->GetComm());
   return glob_dt_est;
}

void C7Operator::ResetVelocityStepEstimate() const
{
   quad_data.dt_est = numeric_limits<double>::infinity();
}

void C7Operator::ComputeDensity(ParGridFunction &rho)
{
   rho.SetSpace(&L2FESpace);

   DenseMatrix Mrho(l2dofs_cnt);
   Vector rhs(l2dofs_cnt), rho_z(l2dofs_cnt);
   Array<int> dofs(l2dofs_cnt);
   DenseMatrixInverse inv(&Mrho);
   MassIntegrator mi(&integ_rule);
   DensityIntegrator di(quad_data);
   di.SetIntRule(&integ_rule);
   for (int i = 0; i < nzones; i++)
   {
      di.AssembleRHSElementVect(*L2FESpace.GetFE(i),
                                *L2FESpace.GetElementTransformation(i), rhs);
      mi.AssembleElementMatrix(*L2FESpace.GetFE(i),
                               *L2FESpace.GetElementTransformation(i), Mrho);
      inv.Factor();
      inv.Mult(rhs, rho_z);
      L2FESpace.GetElementDofs(i, dofs);
      rho.SetSubVector(dofs, rho_z);
   }
}

void C7Operator::PrintTimingData(bool IamRoot, int steps)
{
   double my_rt[5], rt_max[5];
   my_rt[0] = timer.sw_cgH1.RealTime();
   my_rt[1] = timer.sw_cgL2.RealTime();
   my_rt[2] = timer.sw_force.RealTime();
   my_rt[3] = timer.sw_qdata.RealTime();
   my_rt[4] = my_rt[0] + my_rt[2] + my_rt[3];
   MPI_Reduce(my_rt, rt_max, 5, MPI_DOUBLE, MPI_MAX, 0, H1FESpace.GetComm());

   double mydata[2], alldata[2];
   mydata[0] = timer.L2dof_iter;
   mydata[1] = timer.quad_tstep;
   MPI_Reduce(mydata, alldata, 2, MPI_DOUBLE, MPI_SUM, 0, H1FESpace.GetComm());

   if (IamRoot)
   {
      using namespace std;
      cout << endl;
      cout << "CG (H1) total time: " << rt_max[0] << endl;
      cout << "CG (H1) rate (megadofs x cg_iterations / second): "
           << 1e-6 * timer.H1dof_iter / rt_max[0] << endl;
      cout << endl;
      cout << "CG (L2) total time: " << rt_max[1] << endl;
      cout << "CG (L2) rate (megadofs x cg_iterations / second): "
           << 1e-6 * alldata[0] / rt_max[1] << endl;
      cout << endl;
      cout << "Divergences total time: " << rt_max[2] << endl;
      cout << "Divergences rate (megadofs x timesteps / second): "
           << 1e-6 * timer.dof_tstep / rt_max[2] << endl;
      cout << endl;
      cout << "UpdateQuadData total time: " << rt_max[3] << endl;
      cout << "UpdateQuadData rate (megaquads x timesteps / second): "
           << 1e-6 * alldata[1] / rt_max[3] << endl;
      cout << endl;
      cout << "Major kernels total time (seconds): " << rt_max[4] << endl;
      cout << "Major kernels total rate (megadofs x time steps / second): "
           << 1e-6 * H1FESpace.GlobalTrueVSize() * steps / rt_max[4] << endl;
   }
}

void C7Operator::UpdateQuadratureData(double velocity, const Vector &S) const
{
   if (quad_data_is_current) { return; }
   timer.sw_qdata.Start();

   const double N_x = AWBSPhysics->GetVelocityMultiple(); 
   const double N_x_vTmax = AWBSPhysics->GetVelocityScale();
   const double velocity_real = velocity * N_x_vTmax;
   const int nqp = integ_rule.GetNPoints();

   AWBSPhysics->mspei_pcf->SetVelocityReal(velocity_real);
   AWBSPhysics->mspee_pcf->SetVelocityReal(velocity_real);

   ParGridFunction F0, F1;
   Vector* sptr = (Vector*) &S;
   F0.MakeRef(&L2FESpace, *sptr, 0);
   F1.MakeRef(&H1FESpace, *sptr, L2FESpace.GetVSize());

   Vector vector_vals(h1dofs_cnt * dim);
   DenseMatrix Jpi(dim), Jinv(dim), F0stress(dim), F0stressJiT(dim),
               F1stress(dim), F1stressJiT(dim),
               vecvalMat(vector_vals.GetData(), h1dofs_cnt, dim);
   Array<int> L2dofs, H1dofs;

   // Isotropic unit matrix.
   DenseMatrix I; 
   I.Diag(1.0, dim);
   // Isotropic P1 matrix.
   DenseMatrix P1;
   P1.Diag(1.0 / 3.0, dim);

   // Batched computations are needed, because hydrodynamic codes usually
   // involve expensive computations of material properties. Although this
   // miniapp uses simple EOS equations, we still want to represent the batched
   // cycle structure.
   int nzones_batch = 3;
   const int nbatches =  nzones / nzones_batch + 1; // +1 for the remainder.
   int nqp_batch = nqp * nzones_batch;
   double *mspee_b = new double[nqp_batch],
          *mspei_b = new double[nqp_batch],
		  *rho_b = new double[nqp_batch],
          *vTe_b = new double[nqp_batch];
   // Jacobians of reference->physical transformations for all quadrature
   // points in the batch.
   DenseTensor *Jpr_b = new DenseTensor[nqp_batch],
               *AM1_b = new DenseTensor[nqp_batch];
   for (int b = 0; b < nbatches; b++)
   {
      int z_id = b * nzones_batch; // Global index over zones.
      // The last batch might not be full.
      if (z_id == nzones) { break; }
      else if (z_id + nzones_batch > nzones)
      {
         nzones_batch = nzones - z_id;
         nqp_batch    = nqp * nzones_batch;
      }

      for (int z = 0; z < nzones_batch; z++)
      {
         ElementTransformation *T = H1FESpace.GetElementTransformation(z_id);
		 Jpr_b[z].SetSize(dim, dim, nqp);
         AM1_b[z].SetSize(dim, dim, nqp);

         for (int q = 0; q < nqp; q++)
         {
            const IntegrationPoint &ip = integ_rule.IntPoint(q);
            T->SetIntPoint(&ip);
            Jpr_b[z](q) = T->Jacobian(); 

            const double detJ = Jpr_b[z](q).Det();

            const int idx = z * nqp + q;
            //rho_b[idx] = mspei_pcf->GetRho(*T, ip);
            rho_b[idx] = quad_data.rho0DetJ0w(z_id*nqp + q) / 
                         detJ / ip.weight;
            mspei_b[idx] = AWBSPhysics->mspei_pcf->Eval(*T, ip, rho_b[idx]);
            mspee_b[idx] = AWBSPhysics->mspee_pcf->Eval(*T, ip, rho_b[idx]);
            vTe_b[idx] = AWBSPhysics->mspee_pcf->GetvTe(*T, ip); 

            // M1 closure.
            // Matric closure maximizing angular entropy
            // A = 1/3*I + M^2/2*(1 + M^2)*((f1xf1^T)/f1^2 - 1/3*I),
            // where M = |f1|/|f0| must be in (0, 1), "isotropic-freestreaming".
            if (M1_closure)
            {
               double normlim = 1e-32;
               double anisolim = 1e-1;
               //double f0norm = F0.GetValue((*T).ElementNo, ip);
               double f0norm = abs(F0.GetValue((*T).ElementNo, ip));
               //double f0norm = max(normlim, F0.GetValue((*T).ElementNo, ip));
               Vector f1;
               F1.GetVectorValue((*T).ElementNo, ip, f1);
               double f1norm = f1.Norml2();

               if (f0norm < normlim || f1norm < normlim)
               {
                  AM1_b[z](q) = P1;
               }
               //else if (f1norm / f0norm < anisolim)
               //{
               //   AM1_b[z](q) = P1;
               //}
               else
               {
                  double M = min(f1norm / f0norm, 1.0);
                  double Msquare = M * M;
                  double c = Msquare / 2.0 * (1.0 + Msquare);
				  cout << "(f1/f0)^2 / 2.0 * (1.0 + (f1/f0)^2): " << c << endl 
				       << flush;
                  // Shotr version for 1D.
                  AM1_b[z](q) = P1;
                  AM1_b[z](q) *= 1.0 + 2.0 * c;
				  /*
				  DenseMatrix normf1xf1T(dim);
                  normf1xf1T.Diag(1.0, dim);
                  //f1 *= 1.0 / f1norm;
			      // f1 directional matrix.
                  //MultVVt(f1, normf1xf1T);
                  //normf1xf1T = 0.0;
			      //for (int vd = 0; vd < dim; vd++)
                  //{
                  //   normf1xf1T(vd, vd) = f1(vd) * f1(vd);
                  //}
                  // Construct the closure matrix.
			      AM1_b[z](q) = 0.0;
                  AM1_b[z](q).Add(1.0 - c, P1);
                  AM1_b[z](q).Add(c, normf1xf1T);
                  */
				  //cout << "f0norm, f1norm, c:" << f0norm << ", " << f1norm
			      //     << ", " << c << endl << flush;
                  //cout << "normf1xf1T:" << endl << flush;
                  //normf1xf1T.Print();
                  //cout << "AM1:" << endl << flush;
                  //AM1_b[z](q).Print();
               }
            }
         }
         ++z_id;
      }

      z_id -= nzones_batch;
      for (int z = 0; z < nzones_batch; z++)
      {
         ElementTransformation *T = H1FESpace.GetElementTransformation(z_id);
         for (int q = 0; q < nqp; q++)
         {
            const IntegrationPoint &ip = integ_rule.IntPoint(q);
            T->SetIntPoint(&ip);
            //double f0min = 1e-32;
            //double f0 = max(f0min, F0.GetValue((*T).ElementNo, ip));
            //Vector f1;
            //F1.GetVectorValue((*T).ElementNo, ip, f1);

            // Note that the Jacobian was already computed above. We've chosen
            // not to store the Jacobians for all batched quadrature points.
            const DenseMatrix &Jpr = Jpr_b[z](q);
            CalcInverse(Jpr, Jinv);
            const double detJ = Jpr.Det(); 
            const double mspei = mspei_b[z*nqp + q];
            const double mspee = mspee_b[z*nqp + q];
			double rho = rho_b[z*nqp + q];
            double vTe = vTe_b[z*nqp + q];
			// VEF transport closure matrix is either P1 or M1.
			DenseMatrix A1 = P1;
			if (M1_closure) { A1 = AM1_b[z](q); }

            Vector Efield(dim), Bfield(dim), AEfield(dim); 
            AWBSPhysics->Efield_pcf->Eval(Efield, *T, ip);
			AWBSPhysics->Bfield_pcf->Eval(Bfield, *T, ip);

            double Efield_scale;
			AWBSPhysics->Efield_pcf->GetEscales(*T, ip, velocity_real, mspee,
                                                Efield_scale);

			// Matrix projections. 
            A1.Mult(Efield, AEfield); 

            // Stress matrices for f0 and f1 equations.
			F0stress = I;
			// P1 or M1 closure. See the construction above.
            F1stress = A1;

            // Quadrature data for partial assembly of the force operator.
            MultABt(F0stress, Jinv, F0stressJiT);
            F0stressJiT *= integ_rule.IntPoint(q).weight * detJ;
            MultABt(F1stress, Jinv, F1stressJiT);
            F1stressJiT *= integ_rule.IntPoint(q).weight * detJ;
            for (int vd = 0 ; vd < dim; vd++)
            {
               for (int gd = 0; gd < dim; gd++)
               {
                  quad_data.stress1JinvT(vd)(z_id*nqp + q, gd) =
                     F1stressJiT(vd, gd);
                  quad_data.stress0JinvT(vd)(z_id*nqp + q, gd) =
                     F0stressJiT(vd, gd);
               }
               // Extensive vector quadrature data.
               // F0 equation Scattering Efield quadrature 
               // plus compensation of Maxwellization Efield.
               quad_data.Escatter_invrho(z_id*nqp + q, vd) = Efield(vd) / rho;
               quad_data.Escatter_invrho(z_id*nqp + q, vd) += (1.0 
                                                              - Efield_scale) 
                                                              * Efield(vd) 
                                                              / rho;
               // F0 Maxwellization Efield under effect of scaling.
               quad_data.Escaled_invrho(z_id*nqp + q, vd) = Efield_scale * 
                                                            Efield(vd) / rho;
               // F1 full Efield quadrature.
               quad_data.AE_invrho(z_id*nqp + q, vd) = AEfield(vd) / rho;
               // F1 Maxwellization Efield under effect of scaling.
               quad_data.AEscaled_invrho(z_id*nqp + q, vd) = Efield_scale *
                                                             AEfield(vd) / rho;
               // F1 Full Bfield quadrature.
               quad_data.B_invrho(z_id*nqp + q, vd) = Bfield(vd) / rho;
            }

            // Extensive scalar quadrature data.
            quad_data.nu_invrho(z_id*nqp + q) = mspee / rho; //nue/rho; 
            // Scattering on ions and electrons.
			quad_data.nut_invrho(z_id*nqp + q) = (mspei + mspee) / rho;

            // Time step estimate at the point. Here the more relevant length
            // scale is related to the actual mesh deformation; we use the min
            // singular value of the ref->physical Jacobian. In addition, the
            // time step estimate should be aware of the presence of shocks.
            const double h_min =
               Jpr.CalcSingularvalue(dim-1) / (double) H1FESpace.GetOrder(0);

			double dv = h_min * min(mspee, mspei) / N_x_vTmax;;
            quad_data.dt_est = min(quad_data.dt_est, cfl * dv);		
         }
         ++z_id;
      }
   }
   delete [] mspee_b;
   delete [] mspei_b;
   delete [] rho_b;
   delete [] vTe_b;
   delete [] Jpr_b;
   delete [] AM1_b; 
   quad_data_is_current = true;

   timer.sw_qdata.Stop();
   timer.quad_tstep += nzones * nqp;
}

} // namespace nth

} // namespace mfem

#endif // MFEM_USE_MPI

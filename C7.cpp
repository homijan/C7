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
//
//                     __                __
//                    / /   ____  ____  / /_  ____  _____
//                   / /   / __ `/ __ `/ __ \/ __ \/ ___/
//                  / /___/ /_/ / /_/ / / / / /_/ (__  )
//                 /_____/\__,_/\__, /_/ /_/\____/____/
//                             /____/
//
//             High-order Lagrangian Hydrodynamics Miniapp
//
// Laghos(LAGrangian High-Order Solver) is a miniapp that solves the
// time-dependent Euler equation of compressible gas dynamics in a moving
// Lagrangian frame using unstructured high-order finite element spatial
// discretization and explicit high-order time-stepping. Laghos is based on the
// numerical algorithm described in the following article:
//
//    V. Dobrev, Tz. Kolev and R. Rieben, "High-order curvilinear finite element
//    methods for Lagrangian hydrodynamics", SIAM Journal on Scientific
//    Computing, (34) 2012, pp.B606–B641, https://doi.org/10.1137/120864672.
//
// Sample runs:
//    mpirun -np 8 laghos -p 0 -m data/square01_quad.mesh -rs 3 -tf 0.75
//    mpirun -np 8 laghos -p 0 -m data/square01_tri.mesh  -rs 1 -tf 0.75
//    mpirun -np 8 laghos -p 0 -m data/cube01_hex.mesh    -rs 1 -tf 2.0
//    mpirun -np 8 laghos -p 1 -m data/square01_quad.mesh -rs 3 -tf 0.8
//    mpirun -np 8 laghos -p 1 -m data/square01_quad.mesh -rs 0 -tf 0.8 -ok 7 -ot 6
//    mpirun -np 8 laghos -p 1 -m data/cube01_hex.mesh    -rs 2 -tf 0.6
//    mpirun -np 8 laghos -p 2 -m data/segment01.mesh     -rs 5 -tf 0.2
//    mpirun -np 8 laghos -p 3 -m data/rectangle01_quad.mesh -rs 2 -tf 2.5
//    mpirun -np 8 laghos -p 3 -m data/box01_hex.mesh        -rs 1 -tf 2.5
//
// Test problems:
//    p = 0  --> Taylor-Green vortex (smooth problem).
//    p = 1  --> Sedov blast.
//    p = 2  --> 1D Sod shock tube.
//    p = 3  --> Triple point.


#include "laghos_solver.hpp"
#include "c7_solver.hpp"
#include "eos.hpp"
#include "ic.hpp"
#include <memory>
#include <iostream>
#include <fstream>

using namespace std;
using namespace mfem;
using namespace mfem::hydrodynamics;

// Choice for the problem setup.
// int problem; 

void display_banner(ostream & os);

int main(int argc, char *argv[])
{
   // Initialize MPI.
   MPI_Session mpi(argc, argv);
   int myid = mpi.WorldRank();

   // Print the banner.
   if (mpi.Root()) { display_banner(cout); }

   // Parse command-line options.
   const char *mesh_file = "data/square01_quad.mesh";
   int rs_levels = 0;
   int rp_levels = 0;
   int order_v = 2;
   int order_e = 1;
   int ode_solver_type = 4;
   double t_final = 0.5;
   double cfl = 0.5;
   double cg_tol = 1e-8;
   int cg_max_iter = 13;
   int max_tsteps = -1;
   bool p_assembly = true;
   bool visualization = false;
   int vis_steps = 10;
   bool visit = false;
   bool gfprint = false;
   const char *basename = "results/tmp/C7";
   int nth_problem = 5;
   double T_max = 1000.0, T_min = 100.0, rho_max = 10.0, rho_min = 1.0;
   double L = 1.0, T_gradscale = 50.0, rho_gradscale = 50.0;
   double sigma = 8.1027575e17;
   double coulLog = 10.0;
   double Zbar = 0.0; // 4.0;
   // Correct input for a Lorentz force calculation with SH Efield.
   //double EfieldS0 = 1.0;
   //double F0SourceS0 = 1.0;
   // Appropriate value to mimic exactly the SH Efield effect on qH.
   double F0SourceS0 = 0.2857142857142857;
   double EfieldS0 = 0.0;
   double nuS0 = -1.0; // In the case of < 0.0, use the Zbar correction of nu.
   // We expect rho = 1, and so, the ion mass follows.
   double ni = 0.0; //5e19;
   bool M1closure = false;
   // Minimum number of velocity groups.
   double MinimumGroups = 10.0;
   // The point where the detailed kinetic output will be stored at.
   double xpoint = 0.5;
   double c7cfl = 0.25;
   // Number of consistent Efield iterations.
   double dEit_norm_limit = 1e-2;
   int Efield_consistent_iter_max = 100;
   // Efield effect needs to be divided into a directional and an artificial 
   // frictional parts.
   double times_nue = 1.0;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&rs_levels, "-rs", "--refine-serial",
                  "Number of times to refine the mesh uniformly in serial.");
   args.AddOption(&rp_levels, "-rp", "--refine-parallel",
                  "Number of times to refine the mesh uniformly in parallel.");
   args.AddOption(&nth_problem, "-p", "--problem", "Problem setup to use.");
   args.AddOption(&order_v, "-ok", "--order-kinematic",
                  "Order (degree) of the kinematic finite element space.");
   args.AddOption(&order_e, "-ot", "--order-thermo",
                  "Order (degree) of the thermodynamic finite element space.");
   args.AddOption(&ode_solver_type, "-s", "--ode-solver",
                  "ODE solver: 1 - Backward Euler, 2 - SDIRK23 SSP,\n\t" 
				  "            3 - SDIRK 33,       4 - SDIRK 43,\n\t"
                  "            5 - Forward Euler,  6 - RK4.");
   args.AddOption(&t_final, "-tf", "--t-final",
                  "Final time; start time is 0.");
   args.AddOption(&c7cfl, "-cfl", "--cfl", "CFL-condition number.");
   args.AddOption(&cg_tol, "-cgt", "--cg-tol",
                  "Relative CG tolerance (velocity linear solve).");
   args.AddOption(&cg_max_iter, "-cgm", "--cg-max-steps",
                  "Maximum number of CG iterations (velocity linear solve).");
   args.AddOption(&max_tsteps, "-ms", "--max-steps",
                  "Maximum number of steps (negative means no restriction).");
   args.AddOption(&p_assembly, "-pa", "--partial-assembly", "-fa",
                  "--full-assembly",
                  "Activate 1D tensor-based assembly (partial assembly).");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&vis_steps, "-vs", "--visualization-steps",
                  "Visualize every n-th timestep.");
   args.AddOption(&visit, "-visit", "--visit", "-no-visit", "--no-visit",
                  "Enable or disable VisIt visualization.");
   args.AddOption(&gfprint, "-print", "--print", "-no-print", "--no-print",
                  "Enable or disable result output (files in mfem format).");
   args.AddOption(&basename, "-k", "--outputfilename",
                  "Name of the visit dump files");
   args.AddOption(&sigma, "-sigma", "--sigma",
                  "Mean-free-path scaling, i.e. lambda = v^4/sigma/coulLog/ni/Zbar/Zbar.");
   args.AddOption(&coulLog, "-cl", "--coulLog",
                  "Coulomb logarithm used in mean-free-path, i.e. lambda = v^4/sigma/coulLog/ni/Zbar/Zbar.");
   args.AddOption(&Zbar, "-Z", "--Zbar",
                  "Constant ionization used for nu_ee. Used along IGEOS only.");
   args.AddOption(&T_max, "-Tmax", "--Tmax",
                  "Maximum temperature in the step function tanh(x).");
   args.AddOption(&T_min, "-Tmin", "--Tmin",
                  "Minimum temperature in the step function tanh(x).");
   args.AddOption(&rho_max, "-rmax", "--rhomax",
                  "Maximum density in the step function tanh(x).");
   args.AddOption(&rho_min, "-rmin", "--rhomin",
                  "Minimum density in the step function tanh(x).");
   args.AddOption(&L, "-L", "--Length",
                  "Unit interval scale in the function tanh(a*x).");
   args.AddOption(&T_gradscale, "-Tgrad", "--Tgrad",
                  "Temperature gradient scale in the function tanh(a*x).");
   args.AddOption(&rho_gradscale, "-rgrad", "--rhograd",
                  "Density gradient scale in the function tanh(a*x).");
   args.AddOption(&dEit_norm_limit, "-dE", "--dEfielditerationnorm",
                  "Change in the norm related to Electric field to stop iterations.");
   args.AddOption(&Efield_consistent_iter_max, "-Em", "--Emaxit",
                  "Maximum number of the the Electric field iterations.");
   args.AddOption(&EfieldS0, "-E0", "--ES0",
                  "Electric field scaling, i.e. E = S0*E.");
   args.AddOption(&F0SourceS0, "-S0", "--S0",
                  "Electron source scaling (via electron density), i.e. ne = S0*ne.");
   args.AddOption(&nuS0, "-n0", "--n0",
                  "Electron collision frequency scaling, i.e. nu_e = S0*nu_e. If negative Zbar correction is used.");
   args.AddOption(&times_nue, "-xn", "--xnue",
                  "Electric field effect division, maximum friction part as times the nue.");
   args.AddOption(&ni, "-ni", "--ni",
                  "Ion density (conversion as ni = rho/mi).");
   args.AddOption(&M1closure, "-M1", "--M1closure", "-no-M1",
                  "--no-M1closure->P1closure",
                  "Enable or disable M1 VEF closure. If disabled P1 closure applies.");
   args.AddOption(&MinimumGroups, "-minG", "--minimumGroups",
                  "The minimum number of velocity groups.");
   args.AddOption(&xpoint, "-xp", "--xpoint",
                  "An x point within the mesh, where a detailed kinetic profiles as stored.");

   args.Parse();
   if (!args.Good())
   {
      if (mpi.Root()) { args.PrintUsage(cout); }
      return 1;
   }
   if (mpi.Root()) { args.PrintOptions(cout); }

   nth::nth_problem = nth_problem;
   nth::T_max = T_max;
   nth::T_min = T_min;
   nth::rho_max = rho_max;
   nth::rho_min = rho_min;
   nth::L = L;
   nth::T_gradscale = T_gradscale;
   nth::rho_gradscale = rho_gradscale;
   nth::sigma = sigma;
   nth::coulLog = coulLog; // TMP, will be moved to the eos.
   // Read input temperature profile.
   nth::InputProfile inElectronTemp;
   inElectronTemp.SetFile("VFPdata/temperature.dat");
   //inElectronTemp.SetConst(1000.0);
   // Read input density profile. 
   nth::InputProfile inElectronDens; 
   inElectronDens.SetFile("VFPdata/ne.dat");
   if (ni > 0.0) { inElectronDens.SetConst(ni * Zbar); }
   // Read input ionization profile. 
   nth::InputProfile inZbar;
   inZbar.SetFile("VFPdata/zbar.dat");
   if (Zbar > 0.0) { inZbar.SetConst(Zbar); }
   if (nth_problem == 9)
   {
      nth::inElectronTemp = &inElectronTemp;
   }
   if (nth_problem == 10)
   {
      nth::inElectronTemp = &inElectronTemp;
	  nth::inElectronDens = &inElectronDens;
      nth::inZbar = &inZbar;
   }

   // Read the serial mesh from the given mesh file on all processors.
   // Refine the mesh in serial to increase the resolution.
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   const int dim = mesh->Dimension();
   for (int lev = 0; lev < rs_levels; lev++) { mesh->UniformRefinement(); }

   if (p_assembly && dim == 1)
   {
      p_assembly = false;
      if (mpi.Root())
      {
         cout << "Laghos does not support PA in 1D. Switching to FA." << endl;
      }
   }

   // Parallel partitioning of the mesh.
   ParMesh *pmesh = NULL;
   const int num_tasks = mpi.WorldSize();
   const int partitions = floor(pow(num_tasks, 1.0 / dim) + 1e-2);
   int *nxyz = new int[dim];
   int product = 1;
   for (int d = 0; d < dim; d++)
   {
      nxyz[d] = partitions;
      product *= partitions;
   }
   if (product == num_tasks)
   {
      int *partitioning = mesh->CartesianPartitioning(nxyz);
      pmesh = new ParMesh(MPI_COMM_WORLD, *mesh, partitioning);
      delete partitioning;
   }
   else
   {
      if (myid == 0)
      {
         cout << "Non-Cartesian partitioning through METIS will be used.\n";
#ifndef MFEM_USE_METIS
         cout << "MFEM was built without METIS. "
              << "Adjust the number of tasks to use a Cartesian split." << endl;
#endif
      }
#ifndef MFEM_USE_METIS
      return 1;
#endif
      pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   }
   delete [] nxyz;
   delete mesh;

   // Apply the mesh scaling. 
   // In ic.hpp/cpp the original mesh is considered to be a unit (0, 1).
   *(pmesh->GetNodes()) *= L;

   // Refine the mesh further in parallel to increase the resolution.
   for (int lev = 0; lev < rp_levels; lev++) { pmesh->UniformRefinement(); }

   // Define the parallel finite element spaces. We use:
   // - H1 (Gauss-Lobatto, continuous) for position and velocity.
   // - L2 (Bernstein, discontinuous) for specific internal energy.
   // Fundamental! BasisType::Positive destroys C7 calculation with Efield! 
   L2_FECollection L2FEC(order_e, dim, BasisType::GaussLegendre);
   //L2_FECollection L2FEC(order_e, dim, BasisType::GaussLobatto);  
   //L2_FECollection L2FEC(order_e, dim, BasisType::Positive); 
   H1_FECollection H1FEC(order_v, dim);
   ParFiniteElementSpace L2FESpace(pmesh, &L2FEC);
   ParFiniteElementSpace H1FESpace(pmesh, &H1FEC, pmesh->Dimension());

   // ND contains Nedelec "edge-centered" vector finite elements with continuous
   // tangential component.
   ND_FECollection HCurlFEC(order_v, dim);
   ParFiniteElementSpace HCurlFESpace(pmesh, &HCurlFEC);
   ParGridFunction EfieldNedelec_gf(&HCurlFESpace);
   // RT contains Raviart-Thomas "face-centered" vector finite elements with
   // continuous normal component.
   //RT_FECollection HDivFEC(order_v-1, dim);
   //ParFiniteElementSpace  HDivFESpace(pmesh, &HDivFEC);

   // Boundary conditions: all tests use v.n = 0 on the boundary, and we assume
   // that the boundaries are straight.
   Array<int> ess_tdofs;
   {
      Array<int> ess_bdr(pmesh->bdr_attributes.Max()), tdofs1d;
      for (int d = 0; d < pmesh->Dimension(); d++)
      {
         // Attributes 1/2/3 correspond to fixed-x/y/z boundaries, i.e., we must
         // enforce v_x/y/z = 0 for the velocity components.
         ess_bdr = 0; ess_bdr[d] = 1;
         H1FESpace.GetEssentialTrueDofs(ess_bdr, tdofs1d, d);
         ess_tdofs.Append(tdofs1d);
      }
   }

   // Define the explicit ODE solver used for time integration.
   ODESolver *ode_solver = NULL;
   switch (ode_solver_type)
   {
      case 1: ode_solver = new ForwardEulerSolver; break;
      case 2: ode_solver = new RK2Solver; break; 
      case 3: ode_solver = new RK3SSPSolver; break;
      case 4: ode_solver = new RK4Solver; break;
      case 5: ode_solver = new RK2Solver(0.5); break;
      case 6: ode_solver = new RK6Solver; break;
      default:
         if (myid == 0)
         {
            cout << "Unknown ODE solver type: " << ode_solver_type << '\n';
         }
         delete pmesh;
         MPI_Finalize();
         return 3;
   }

   HYPRE_Int glob_size_l2 = L2FESpace.GlobalTrueVSize();
   HYPRE_Int glob_size_h1 = H1FESpace.GlobalTrueVSize();

   if (mpi.Root())
   {
      cout << "Number of kinematic (position, velocity) dofs: "
           << glob_size_h1 << endl;
      cout << "Number of specific internal energy dofs: "
           << glob_size_l2 << endl;
   }

   int Vsize_l2 = L2FESpace.GetVSize();
   int Vsize_h1 = H1FESpace.GetVSize();

   // The monolithic BlockVector stores unknown fields as:
   // - 0 -> position
   // - 1 -> velocity
   // - 2 -> specific internal energy

   Array<int> true_offset(4);
   true_offset[0] = 0;
   true_offset[1] = true_offset[0] + Vsize_h1;
   true_offset[2] = true_offset[1] + Vsize_h1;
   true_offset[3] = true_offset[2] + Vsize_l2;
   BlockVector S(true_offset);

   // Define GridFunction objects for the position, velocity and specific
   // internal energy.  There is no function for the density, as we can always
   // compute the density values given the current mesh position, using the
   // property of pointwise mass conservation.
   ParGridFunction x_gf, v_gf, e_gf;
   x_gf.MakeRef(&H1FESpace, S, true_offset[0]);
   v_gf.MakeRef(&H1FESpace, S, true_offset[1]);
   e_gf.MakeRef(&L2FESpace, S, true_offset[2]);

   // Initialize x_gf using the starting mesh coordinates. This also links the
   // mesh positions to the values in x_gf.
   pmesh->SetNodalGridFunction(&x_gf);

   // Initialize the velocity.
   VectorFunctionCoefficient v_coeff(pmesh->Dimension(), nth::v0);
   v_gf.ProjectCoefficient(v_coeff);

   // Initialize density and specific internal energy values. We interpolate in
   // a non-positive basis to get the correct values at the dofs.  Then we do an
   // L2 projection to the positive basis in which we actually compute. The goal
   // is to get a high-order representation of the initial condition. Note that
   // this density is a temporary function and it will not be updated during the
   // time evolution.
   ParGridFunction rho_gf(&L2FESpace);
   FunctionCoefficient rho_coeff(nth::rho0);
   L2_FECollection l2_fec(order_e, pmesh->Dimension());
   ParFiniteElementSpace l2_fes(pmesh, &l2_fec);
   ParGridFunction l2_rho(&l2_fes), l2_e(&l2_fes);
   l2_rho.ProjectCoefficient(rho_coeff);
   rho_gf.ProjectGridFunction(l2_rho);
   if (nth::nth_problem == 1)
   {
      // For the Sedov test, we use a delta function at the origin.
      DeltaCoefficient e_coeff(0, 0, 0.25);
      l2_e.ProjectCoefficient(e_coeff);
	  // Set min temperature to be nonzero.
      //ParGridFunction cnst(&l2_fes);
	  //cnst = 0.0025;
	  //l2_e += cnst;
   }
   else
   {
      FunctionCoefficient e_coeff(nth::e0);
      l2_e.ProjectCoefficient(e_coeff);
   }
   e_gf.ProjectGridFunction(l2_e);

   // Space-dependent ideal gas coefficient over the Lagrangian mesh.
   FunctionCoefficient gamma_cf = (nth::gamma);
   Coefficient *material_pcf = &gamma_cf;

   // Additional details, depending on the problem.
   int source = 0; bool visc;
   switch (nth::nth_problem)
   {
      case 0: if (pmesh->Dimension() == 2) { source = 1; }
         visc = false; break;
      case 1: visc = true; break;
      case 2: visc = true; break;
      case 3: visc = true; break;
      case 4: visc = true; break;
      case 5: visc = true; break;
      case 6: visc = true; break;
      case 7: visc = true; break;
      case 8: visc = true; break;
      case 9: visc = true; break;
      case 10: visc = true; break;
      default: MFEM_ABORT("Wrong problem specification!");
   }

   LagrangianHydroOperator oper(S.Size(), H1FESpace, L2FESpace,
                                ess_tdofs, rho_gf, source, cfl, material_pcf,
                                visc, p_assembly, cg_tol, cg_max_iter);

///////////////////////////////////////////////////////////////
///// C7 nonlocal solver //////////////////////////////////////
///////////////////////////////////////////////////////////////
   // Efield and Bfield grid functions.
   //ParGridFunction Efield_gf(&HCurlFESpace);
   ParGridFunction Efield_gf(&H1FESpace);
   //ParGridFunction Bfield_gf(&HDivFESpace);

   // The monolithic BlockVector stores unknown fields as:
   // - 0 -> isotropic F0 (energy density)
   // - 1 -> anisotropic F1 (flux density)
   Array<int> c7true_offset(8);
   c7true_offset[0] = 0;                           // F0
   c7true_offset[1] = c7true_offset[0] + Vsize_l2; // F1
   c7true_offset[2] = c7true_offset[1] + Vsize_h1; // qH
   c7true_offset[3] = c7true_offset[2] + Vsize_h1; // j
   c7true_offset[4] = c7true_offset[3] + Vsize_h1; // a0
   c7true_offset[5] = c7true_offset[4] + Vsize_l2; // b0
   c7true_offset[6] = c7true_offset[5] + Vsize_h1; // b1
   c7true_offset[7] = c7true_offset[6] + Vsize_h1;
   BlockVector c7F(c7true_offset);

   // Define GridFunction objects for the zero and first moments of
   // the~electron distribution function.
   ParGridFunction F0_gf, F1_gf, hflux_gf, jC_gf, a0_gf, b0_gf, b1_gf;
   F0_gf.MakeRef(&L2FESpace, c7F, c7true_offset[0]);
   F1_gf.MakeRef(&H1FESpace, c7F, c7true_offset[1]);
   hflux_gf.MakeRef(&H1FESpace, c7F, c7true_offset[2]);
   jC_gf.MakeRef(&H1FESpace, c7F, c7true_offset[3]);
   a0_gf.MakeRef(&L2FESpace, c7F, c7true_offset[4]);
   b0_gf.MakeRef(&H1FESpace, c7F, c7true_offset[5]);
   b1_gf.MakeRef(&H1FESpace, c7F, c7true_offset[6]);

   // Define hydrodynamics related coefficients as mean stopping power and
   // source function depending on plasma temperature and density. 
   const double kB = 1.6022e-12, me = 9.1094e-28, qe = 4.8032e-10, 
                pi = 3.14159265359; 
   const double mi = 1.0 / ni; // Expecting rho = 1.0.
   //const double mi = Zbar / ne;
   //// Define the ideal gas equation of state.
   //nth::IGEOS eosig(me, kB);
   //// Use a constant ionization provided by IG eos.
   //eosig.SetZbar(Zbar);
   //// Use a homogeneous ion mass used within IG eos.
   //eosig.SetIonMass(mi);

   // Define an external equation of state. In principle avoids hydro.
   nth::ExtDataEOS eosext(me, kB);
   // Use a constant electron density.
   eosext.SetConstNe(ni * Zbar);
   // Use a constant ionization.
   eosext.SetConstZbar(Zbar);

   //double Te_test = 1e3;
   //cout << "eosext x(T): " << eosext.GetX(0.0, Te_test) << endl;

   nth::EOS *eos = &eosext; //&eosig;

   // Prepare C6 physics.
   nth::ClassicalMeanStoppingPower mspei_cf(rho_gf, e_gf, v_gf, material_pcf,
                                            eos);
   nth::ClassicalAWBSMeanStoppingPower mspee_cf(rho_gf, e_gf, v_gf, 
                                                material_pcf, eos); 
   if (nuS0 > 0.0) { mspee_cf.SetCorrAWBS(nuS0); }
   
   nth::NTHvHydroCoefficient *mspei_pcf = &mspei_cf;
   nth::NTHvHydroCoefficient *mspee_pcf = &mspee_cf;
   nth::AWBSdfMdv dfMdv_cf(rho_gf, e_gf, v_gf, material_pcf, eos);
   nth::NTHvHydroCoefficient *fM_pcf = dfMdv_cf.GetMaxBo();
   nth::NTHvHydroCoefficient *dfMdv_pcf = &dfMdv_cf;

   // Create original Lorentz zero current Efield and set input scales 
   // of the electron source and the Efield.
   nth::LorentzEfield LorentzEfield_cf(pmesh->Dimension(), rho_gf, e_gf, v_gf, 
                                       material_pcf, eos);  
   dfMdv_cf.SetScale0(F0SourceS0);
   LorentzEfield_cf.SetScale0(EfieldS0);
   // Represent Lorentz Efield as a VectorCoefficient.
   VectorCoefficient *LorentzEfield_pvcf = &LorentzEfield_cf;
   VectorCoefficient &Efield_vcf = *LorentzEfield_pvcf;
   // Estimate the Efield by projecting Lorentz to the Efield grid function.
   //Efield_gf.ProjectCoefficient(Efield_vcf);
   Efield_gf = 0.0;

   // Represent Efield by a vector coefficient.
   nth::EfieldCoefficient Efield_gfcf(&Efield_gf);
   Efield_gfcf.SetTimesNue(times_nue);
   nth::EfieldCoefficient *Efield_pcf = &Efield_gfcf;

   // TMP Prepare a fake Bfield vector coefficient.
   Vector vZero(pmesh->Dimension());
   vZero = 0.0;
   VectorConstantCoefficient ZeroBfield_cf(vZero);
   VectorCoefficient *Bfield_pcf = &ZeroBfield_cf;


   // NEW 
   nth::OhmCurrentCoefficient OhmCurrent_cf(pmesh->Dimension(), &a0_gf, &b0_gf,
                                            &b1_gf, Efield_pcf, Bfield_pcf);
   nth::OhmEfieldCoefficient OhmEfield_cf(pmesh->Dimension(), &jC_gf, &a0_gf, 
                                          &b0_gf, &b1_gf, Bfield_pcf);
   // Pure zero current condition. Alpha gives a portion of previous current.
   OhmEfield_cf.SetAlpha(0.0);

   // This object represents physics in AWBS Boltzmann transport model.
   nth::AWBSMasterOfPhysics AWBSPhysics(pmesh->Dimension(), mspei_pcf, 
                                        mspee_pcf, fM_pcf, dfMdv_pcf,
                                        Efield_pcf, Bfield_pcf, eos);

   // Diagnostics.
   // Create the Efield correction coefficient.
   nth::CorrEfieldCoefficient corrEfield_cf(pmesh->Dimension(), mspei_pcf, 
                                            mspee_pcf); 
   corrEfield_cf.SetEfield(Efield_pcf);
   Coefficient * corrEfield_pcf = &corrEfield_cf;
   GridFunction corrEfield_gf(&L2FESpace);

   // Static coefficient defined in c7_solver.hpp.
   //double c7cfl = 0.25;
   //double c7cfl = 0.005;
   vis_steps = 1000000000;
   // ALWAYS calculate on v in (0, 1)
   double vmax = 1.0;
   double vTmultiple = 6.0; //7.0;
   // well, not really, since the lowest v = 0 is singular, so
   //double vmin = 0.001 * vmax;
   double vmin_multiple = 0.01;
   double vmin = vmin_multiple * vmax;
   //double vmin = 0.07 * vmax;
   // and provide some maximum dv step.
   double dvmax = (vmax - vmin) / MinimumGroups;
   //double dvmax = vmax*0.0005;
   bool nonlocal_test = false;
   if (nonlocal_test)
   {
      vTmultiple = 6.0;
      vmin = 0.01 * vmax;
      //vmin = 3.5 * vmax / vTmultiple; // Minimum 3.5*vTh
      dvmax = vmax*0.1;
      if (pmesh->Dimension() == 1)
      { 
         nth::sigma = 2e3;
         vis_steps = 10000;
         c7cfl = 0.5;
      }
      else if (pmesh->Dimension() == 2)
      { 
	     nth::sigma = 1e5; //2e1;
         vis_steps = 10000;
         c7cfl = 1.0;
      }
      else if (pmesh->Dimension() == 3)
      {
         nth::sigma = 5e7;
         vis_steps = 10000;
         c7cfl = 0.5;
      }
   }

   oper.ComputeDensity(rho_gf); 
   double loc_Tmax = e_gf.Max(), glob_Tmax;
   MPI_Allreduce(&loc_Tmax, &glob_Tmax, 1, MPI_DOUBLE, MPI_MAX,
                 pmesh->GetComm());
   AWBSPhysics.SetVelocityScale(vTmultiple, glob_Tmax);
   // Check, must return exactly vTmultiple.
   double N_x_vTmax = AWBSPhysics.GetVelocityScale();

   // Initialize the C7-AWBS operator
   nth::C7Operator c7oper(c7F.Size(), H1FESpace, L2FESpace, ess_tdofs, rho_gf, 
                          c7cfl, &AWBSPhysics, x_gf, e_gf, cg_tol, cg_max_iter);
   // Turn on M1closure.
   if (M1closure) { c7oper.SetM1closure(); }

   // Define the explicit/implicit ODE solver used for velocity integration.
   ODESolver *c7ode_solver = NULL; 

   switch (ode_solver_type)
   {
      // L-stable
      case 1: c7ode_solver = new BackwardEulerSolver; break;
      case 2: c7ode_solver = new SDIRK23Solver(2); break;
      case 3: c7ode_solver = new SDIRK33Solver; break;
      // A-stable
      //c7ode_solver = new ImplicitMidpointSolver;
      //c7ode_solver = new SDIRK23Solver;
      case 4: c7ode_solver = new SDIRK34Solver; break;
      // Explicit
	  case 5: c7ode_solver = new ForwardEulerSolver;
      //c7ode_solver = new RK2Solver(0.5);
      case 6: c7ode_solver = new RK4Solver; break;
      //c7ode_solver = new RK6Solver;
      default:
         if (myid == 0)
         {
            cout << "Unknown ODE solver type: " << ode_solver_type << '\n';
         }
         delete pmesh;
         MPI_Finalize();
         return 3;
   }
   
   c7ode_solver->Init(c7oper);
/*   
   c7oper.ResetVelocityStepEstimate();
   c7oper.ResetQuadratureData();
   c7oper.SetTime(vmax);
   //double dvmax = vmax*0.1;
   double dvmin = min(dvmax, c7oper.GetVelocityStepEstimate(c7F));
   F0_gf = 0.0; F1_gf = 0.0;
   int c7ti = 0;
   double v = vmax;
   double dv = -dvmin;
   //EfieldNedelec_gf.ProjectCoefficient(Efield_cf);
   hflux_gf = 0.0;
   jC_gf = 0.0;
   while (abs(dv) >= abs(dvmin))
   {
      c7ti++;
      c7ode_solver->Step(c7F, v, dv);

      // Perform the integration over velocity space.
      intf0_gf.Add(pow(N_x_vTmax*v, 2.0) * N_x_vTmax*abs(dv), F0_gf);
      //jC_gf.Add(pow(N_x_vTmax*v, 3.0) * N_x_vTmax*abs(dv), F1_gf);
      //hflux_gf.Add(me / 2.0 * pow(N_x_vTmax*v, 5.0) * N_x_vTmax*abs(dv), 
	  //             F1_gf);

      double loc_minF0 = F0_gf.Min(), glob_minF0;
      MPI_Allreduce(&loc_minF0, &glob_minF0, 1, MPI_DOUBLE, MPI_MIN,
                       pmesh->GetComm());
      double loc_maxF0 = F0_gf.Max(), glob_maxF0;
      MPI_Allreduce(&loc_maxF0, &glob_maxF0, 1, MPI_DOUBLE, MPI_MAX,
                       pmesh->GetComm());
      double loc_minF1 = F1_gf.Min(), glob_minF1;
      MPI_Allreduce(&loc_minF1, &glob_minF1, 1, MPI_DOUBLE, MPI_MIN,
                       pmesh->GetComm());
      double loc_maxF1 = F1_gf.Max(), glob_maxF1;
      MPI_Allreduce(&loc_maxF1, &glob_maxF1, 1, MPI_DOUBLE, MPI_MAX,
                       pmesh->GetComm());

      c7oper.ResetVelocityStepEstimate();
      c7oper.ResetQuadratureData();
      c7oper.SetTime(v);
      dv = - min(dvmax, c7oper.GetVelocityStepEstimate(c7F));
      if (v + dv < vmin) { dv = vmin - v; }

      if (mpi.Root())
      {
         cout << fixed;
         cout << "group " << setw(5) << c7ti
                 << ",\tv = " << setw(5) << setprecision(4) << v
                 << ",\tdv = " << setw(5) << setprecision(8) << dv << endl
                 << "[min(f0), max(f0)] = [" << setprecision(17)
                 << glob_minF0 << ",\t" << glob_maxF0 << "]" << endl
                 << "[min(f1), max(f1)] = [" << setprecision(17)
                 << glob_minF1 << ",\t" << glob_maxF1 << "]"
                 << endl;
      }
   }
*/
///////////////////////////////////////////////////////////////
///// C7 nonlocal solver //////////////////////////////////////
///////////////////////////////////////////////////////////////

   socketstream vis_rho, vis_v, vis_e, vis_f0, vis_j, vis_Efield, vis_Kn, 
                vis_hflux;
   char vishost[] = "localhost";
   int  visport   = 19916;

   if (visualization || visit) { oper.ComputeDensity(rho_gf); }

   if (visualization)
   {
      // Make sure all MPI ranks have sent their 'v' solution before initiating
      // another set of GLVis connections (one from each rank):
      MPI_Barrier(pmesh->GetComm());

      vis_rho.precision(8);
      vis_v.precision(8);
      vis_e.precision(8);

      vis_f0.precision(8);
      vis_j.precision(8);
      vis_Efield.precision(8);
	  vis_Kn.precision(8);
      vis_hflux.precision(8);

      int Wx = 0, Wy = 0; // window position
      const int Ww = 350, Wh = 350; // window size
      int offx = Ww+10; // window offsets

      VisualizeField(vis_rho, vishost, visport, rho_gf,
                     "Density", Wx, Wy, Ww, Wh);
      //Wx += offx;
      //VisualizeField(vis_v, vishost, visport, v_gf,
      //               "Velocity", Wx, Wy, Ww, Wh);
      Wx += offx;
      VisualizeField(vis_j, vishost, visport, jC_gf,
                     "Current", Wx, Wy, Ww, Wh);
      Wx += offx;
      VisualizeField(vis_e, vishost, visport, e_gf,
                     "T", Wx, Wy, Ww, Wh);

      Wx = 0;
      Wy +=offx;
      VisualizeField(vis_Efield, vishost, visport, Efield_gf, "Efield", 
                     Wx, Wy, Ww, Wh);
	  //VisualizeField(vis_Kn, vishost, visport, Kn_gf, "Kn", 
      //               Wx, Wy, Ww, Wh);
      Wx += offx;
      VisualizeField(vis_hflux, vishost, visport, hflux_gf,
                     "Heat flux", Wx, Wy, Ww, Wh);
   }

   // Save data for VisIt visualization
   VisItDataCollection visit_dc(basename, pmesh);
   if (visit)
   {
      visit_dc.RegisterField("Density",  &rho_gf);
      visit_dc.RegisterField("Velocity", &v_gf);
      visit_dc.RegisterField("Specific Internal Energy", &e_gf);
      visit_dc.SetCycle(0);
      visit_dc.SetTime(0.0);
      visit_dc.Save();
   }

   // Perform time-integration (looping over the time iterations, ti, with a
   // time-step dt). The object oper is of type LagrangianHydroOperator that
   // defines the Mult() method that used by the time integrators.
   ode_solver->Init(oper);
   oper.ResetTimeStepEstimate();
   double t = 0.0, dt = oper.GetTimeStepEstimate(S), t_old;
   bool last_step = false;
   int steps = 0;
   BlockVector S_old(S);
   for (int ti = 1; !last_step; ti++)
   {
      if (t + dt >= t_final)
      {
         dt = t_final - t;
         last_step = true;
      }
      if (steps == max_tsteps) { last_step = true; }

      S_old = S;
      t_old = t;
      oper.ResetTimeStepEstimate();

      // S is the vector of dofs, t is the current time, and dt is the time step
      // to advance.
      ode_solver->Step(S, t, dt);
      steps++;

      // Adaptive time step control.
      const double dt_est = oper.GetTimeStepEstimate(S);
      if (dt_est < dt)
      {
         // Repeat (solve again) with a decreased time step - decrease of the
         // time estimate suggests appearance of oscillations.
         dt *= 0.85;
         if (dt < numeric_limits<double>::epsilon())
         { MFEM_ABORT("The time step crashed!"); }
         t = t_old;
         S = S_old;
         oper.ResetQuadratureData();
         if (mpi.Root()) { cout << "Repeating step " << ti << endl; }
         ti--; continue;
      }
      else if (dt_est > 1.25 * dt) { dt *= 1.02; }

      // Make sure that the mesh corresponds to the new solution state.
      pmesh->NewNodes(x_gf, false);

      if (last_step || (ti % vis_steps) == 0)
      {
         double loc_norm = e_gf * e_gf, tot_norm;
         MPI_Allreduce(&loc_norm, &tot_norm, 1, MPI_DOUBLE, MPI_SUM,
                       pmesh->GetComm());
         if (mpi.Root())
         {
            cout << fixed;
            cout << "step " << setw(5) << ti
                 << ",\tt = " << setw(5) << setprecision(4) << t
                 << ",\tdt = " << setw(5) << setprecision(6) << dt
                 << ",\t|e| = " << setprecision(10)
                 << sqrt(tot_norm) << endl;
         }

         // Make sure all ranks have sent their 'v' solution before initiating
         // another set of GLVis connections (one from each rank):
         MPI_Barrier(pmesh->GetComm());

///////////////////////////////////////////////////////////////
///// C7 nonlocal solver //////////////////////////////////////
///////////////////////////////////////////////////////////////
         oper.ComputeDensity(rho_gf); 
         double loc_Tmax = e_gf.Max(), glob_Tmax;
         MPI_Allreduce(&loc_Tmax, &glob_Tmax, 1, MPI_DOUBLE, MPI_MAX,
                       pmesh->GetComm());
         AWBSPhysics.SetVelocityScale(vTmultiple, glob_Tmax);         
		 // Check, must return exactly vTmultiple.
		 N_x_vTmax = AWBSPhysics.GetVelocityScale();

         // Point value structures for storing the distribution function.
         // For a given point.
		 //int mpirank_point;
		 int cell_point = 0;
         IntegrationPoint ip_point;
         ip_point.Set3(0.5, 0.5, 0.5);
         double f0_point;
         Vector f1_point; 
         vector<double> v_point, f0_v_point, f1x_v_point, f0v2_v_point,
                        mehalff1xv5_v_point, mehalff0v5_v_point; 
         bool right_proc_point = false;
		 // For a point of the nonlocal heat flux maximum.
		 int mpirank_pointmax;
		 int cell_pointmax = 0;
         IntegrationPoint ip_pointmax;
         ip_pointmax.Set3(0.5, 0.5, 0.5);
         double xpointmax;
		 double f0_pointmax;
         Vector f1_pointmax; 
         vector<double> v_pointmax, f0_v_pointmax, f1x_v_pointmax, 
                        f0v2_v_pointmax, mehalff1xv5_v_pointmax, 
                        mehalff0v5_v_pointmax;
         bool right_proc_pointmax = false;
		 // Find an element where the xpoint belongs and find its ip.
         IntegrationPoint ip_min, ip_max;
         ip_min.Set3(0.0, 0.0, 0.0);
         ip_max.Set3(1.0, 1.0, 1.0); 
         int elNo = 0; 
		 while (!right_proc_point & elNo < x_gf.FESpace()->GetNE() - 1)
         {
            double x_min = x_gf.GetValue(elNo, ip_min);
            double x_max = x_gf.GetValue(elNo, ip_max);
            if (xpoint >= x_min & xpoint <= x_max) 
            { 
               // We are on the right processor.
			   right_proc_point = true;
			   // The correct cell number
			   cell_point = elNo;
			   // Find an ip corresponding to xpoint.
               double tol = 1e-3;
               double L = x_max - x_min, dx, xc, sc, sL = 0.0, sR = 1.0;
			   cout << "xpoint: " << xpoint << endl << flush;
			   do
               {
                  sc = (sL + sR) / 2.0;
				  ip_point.Set3(sc, sc, sc);
				  xc = x_gf.GetValue(cell_point, ip_point);
                  if (xc > xpoint) { sR = sc; }
                  else { sL = sc; }
				  dx = abs(xc - xpoint);
				  //cout << "xpoint, xc: " << xpoint << ", " << xc << ", " 
                  //     << sL << ", " << sR << endl << flush; 
			   } while (dx/L > tol);
               // The maximum is yet unknown, so we use xpoint values.
		       //cell_pointmax = cell_point;
		       //ip_pointmax = ip_point;
               //right_proc_pointmax = true;
            }
            //cout << "right_proc_point, elNo, x_min, x_max: " 
            //     << right_proc_point << ", " << elNo << ", " << x_min << ", " 
            //     << x_max << endl << flush;	
			elNo++;
         }

         // Starting value of the E field norm.
		 double loc_Eit_norm = Efield_gf.Norml2();
         double glob_old_Eit_norm = 1.0;
         //MPI_Allreduce(&loc_Eit_norm, &glob_old_Eit_norm, 1, 
         //              MPI_DOUBLE, MPI_SUM, pmesh->GetComm());
         double dEit_norm = 1e64;
         bool converging = true;
		 int Eit = 0;
  
         double sign = -1.0;
         double v_lim;
		 if (sign < 0.0) { v_lim = vmin; }
         else v_lim = vmax;

         while (Eit < Efield_consistent_iter_max && converging)
		 { 
		    // Actual integration of C7Operator.
            c7oper.ResetVelocityStepEstimate();
            c7oper.ResetQuadratureData();
            c7oper.SetTime(vmax);
            double dvmin = min(dvmax, c7oper.GetVelocityStepEstimate(c7F));
            F0_gf = 0.0; //1e-2; 
		    F1_gf = 0.0;
            int c7ti = 0;
            double v; 
			if (sign < 0.0) { v = vmax; }
			else v = vmin;
            double dv = sign * dvmin;
            hflux_gf = 0.0;
            jC_gf = 0.0;
		    a0_gf = 0.0;
		    b0_gf = 0.0;
		    b1_gf = 0.0;	
            // Point distribution function output restart.
            v_point.resize(0); f0_v_point.resize(0); f1x_v_point.resize(0); 
            f0v2_v_point.resize(0); mehalff1xv5_v_point.resize(0); 
            mehalff0v5_v_point.resize(0);
            // Heat flux max point distribution function output restart. 
            v_pointmax.resize(0); f0_v_pointmax.resize(0); 
            f1x_v_pointmax.resize(0); f0v2_v_pointmax.resize(0); 
            mehalff1xv5_v_pointmax.resize(0); mehalff0v5_v_pointmax.resize(0);
			while (sign*(v - v_lim) < 0.0)
		    {
               c7ti++;
               c7ode_solver->Step(c7F, v, dv);
            
               // Store the distribution function at a given point.
               if (right_proc_point)
			   {
                  //cout << "cell_point: " << cell_point << endl << flush;
			      f0_point = F0_gf.GetValue(cell_point, ip_point);
                  F1_gf.GetVectorValue(cell_point, ip_point, f1_point);
                  v_point.push_back(N_x_vTmax * v);
			      f0_v_point.push_back(f0_point);
			      // TODO extend to more dimensions.
			      f1x_v_point.push_back(f1_point(0));
			      f0v2_v_point.push_back(f0_point * pow(N_x_vTmax*v, 2.0));
			      // TODO extend to more dimensions.
			      mehalff1xv5_v_point.push_back(0.5 * me * f1_point(0) * 
                                                pow(N_x_vTmax*v, 5.0));
			      mehalff0v5_v_point.push_back(0.5 * me * f0_point * 
                                                pow(N_x_vTmax*v, 5.0));
               }
               // Store the distribution function at the max heat flux point.
               if (right_proc_pointmax)
			   {
                  //cout << "cell_point: " << cell_point << endl << flush;
			      f0_point = F0_gf.GetValue(cell_pointmax, ip_pointmax);
                  F1_gf.GetVectorValue(cell_pointmax, ip_pointmax, f1_point);
                  v_pointmax.push_back(N_x_vTmax * v);
			      f0_v_pointmax.push_back(f0_point);
			      // TODO extend to more dimensions.
			      f1x_v_pointmax.push_back(f1_point(0));
			      f0v2_v_pointmax.push_back(f0_point * pow(N_x_vTmax*v, 2.0));
			      // TODO extend to more dimensions.
			      mehalff1xv5_v_pointmax.push_back(0.5 * me * f1_point(0) * 
                                                   pow(N_x_vTmax*v, 5.0));
			      mehalff0v5_v_pointmax.push_back(0.5 * me * f0_point * 
                                                  pow(N_x_vTmax*v, 5.0));
               }


               c7oper.ResetVelocityStepEstimate();
               c7oper.ResetQuadratureData();
               c7oper.SetTime(v);
               dv = sign * min(dvmax, c7oper.GetVelocityStepEstimate(c7F));
               if (sign * (v + dv - v_lim) > 0.0) { dv = v_lim - v; }

			   double loc_minF0 = F0_gf.Min(), glob_minF0;
               MPI_Allreduce(&loc_minF0, &glob_minF0, 1, MPI_DOUBLE, MPI_MIN,
                             pmesh->GetComm());
               double loc_maxF0 = F0_gf.Max(), glob_maxF0;
               MPI_Allreduce(&loc_maxF0, &glob_maxF0, 1, MPI_DOUBLE, MPI_MAX,
                             pmesh->GetComm());
               double loc_minF1 = F1_gf.Min(), glob_minF1;
               MPI_Allreduce(&loc_minF1, &glob_minF1, 1, MPI_DOUBLE, MPI_MIN,
                             pmesh->GetComm());
               double loc_maxF1 = F1_gf.Max(), glob_maxF1;
               MPI_Allreduce(&loc_maxF1, &glob_maxF1, 1, MPI_DOUBLE, MPI_MAX,
                             pmesh->GetComm());

               if (mpi.Root())
               {
                  cout << fixed;
                  cout << "group " << setw(5) << c7ti
                  << ",\tv = " << setw(5) << setprecision(4) << v
                  << ",\tdv = " << setw(5) << setprecision(8) << dv << endl
                  << "[min(f0), max(f0)] = [" << setprecision(17)
                  << glob_minF0 << ",\t" << glob_maxF0 << "]" << endl
                  << "[min(f1), max(f1)] = [" << setprecision(17)
                  << glob_minF1 << ",\t" << glob_maxF1 << "]"
                  << endl;
               }
            }
            // After the integration:
			// Since a deceleration (int_vmax->vmin) run has been performed, 
			// the integrated quantities needs to revert sign.
            hflux_gf *= sign;
            jC_gf *= sign;
            // Microscopic Ohm's law / Ampere's law components too.
			a0_gf *= sign;
            b0_gf *= sign;
            b1_gf *= sign;

            // Consistent Efield (zero current) computation.
            Efield_gf.ProjectCoefficient(OhmEfield_cf);

            // TMP testing of jC current convergence.
            //jC_gf.ProjectCoefficient(OhmCurrent_cf);
            // Gather the current norm.
            double loc_jC_norm = jC_gf.Norml2();
			double glob_jC_norm;
            MPI_Allreduce(&loc_jC_norm, &glob_jC_norm, 1, 
                          MPI_DOUBLE, MPI_SUM, pmesh->GetComm());
            // Gather the heat flux norm.
			double loc_hflux_norm = hflux_gf.Norml2();
			double glob_hflux_norm;
            MPI_Allreduce(&loc_hflux_norm, &glob_hflux_norm, 1, 
                          MPI_DOUBLE, MPI_SUM, pmesh->GetComm());
            // Use the current norm as converging criteria.
			double glob_new_Eit_norm = glob_jC_norm;
            // Use the heat flux norm as converging criteria.
			//double glob_new_Eit_norm = glob_hflux_norm;
            double dEit_norm_new = abs(glob_old_Eit_norm - glob_new_Eit_norm) 
                                  / glob_old_Eit_norm;
            // Stop iteration if reached convergence or not converging.
            if (dEit_norm_new < dEit_norm_limit)
            //if (dEit_norm_new > dEit_norm || dEit_norm_new < dEit_norm_limit)
            { 
               converging = false; 
            }
            dEit_norm = dEit_norm_new;
            if (mpi.Root())
            {
               cout << "Eit, jC_norm, dhflux_norm(|j0-j1|/|j0|): " << Eit 
                    << ", " << std::scientific << setprecision(4) 
                    << glob_jC_norm << ", " << dEit_norm << endl;
            }
            glob_old_Eit_norm = glob_new_Eit_norm;

            // Maximum flux point look up.
            // Construct a double-int memory structure.
			struct { 
               double val; 
               int rank; 
            } loc_hflux_doubleint, glob_hflux_doubleint;
            // Store values.
			loc_hflux_doubleint.val = hflux_gf.Max();
			loc_hflux_doubleint.rank = myid;
			// Reduce the correct mpirank.
			MPI_Allreduce(&loc_hflux_doubleint, &glob_hflux_doubleint, 1, 
                          MPI_DOUBLE_INT, MPI_MAXLOC, pmesh->GetComm());
			// Go over all cells and find pointmax.
            mpirank_pointmax = glob_hflux_doubleint.rank;
			if(myid == mpirank_pointmax)
            {
               double hflux_max = 0.0;
			   IntegrationPoint loc_ip_point;
			   for (int elNo = 0; elNo < x_gf.FESpace()->GetNE(); elNo++)
               {
                  int N = 10;
                  for (int p = 0; p < N; p++)
                  {
                     double sc = p * 1./(N-1);
                     loc_ip_point.Set3(sc, sc, sc);
                     double hflux = hflux_gf.GetValue(elNo, loc_ip_point);
					 if (hflux > hflux_max)
                     {
                        hflux_max = hflux;
                        // Values used to evaluate and store grid functions.
                        cell_pointmax = elNo; 
                        ip_pointmax = loc_ip_point;
                        // Evaluate x coordinate.
						xpointmax = x_gf.GetValue(elNo, ip_pointmax);
                     }
                  }
               }
               // Allow for distribution function storing.
			   right_proc_pointmax = true;

               cout << "xpointmax, hflux_max: " << xpointmax << ", " 
                    << 1e-7 * 0.5 * me * hflux_max << endl;
            }
            else { right_proc_pointmax = false; }

            // Treat internally integrated entities.
            // Heat flux.
            hflux_gf *= me / 2.0; // c7oper does not use me.
            // Current.
            jC_gf *= qe; // c7oper does not use qe.
            // End of this loop.
		    Eit++;
         } // Efield loop.

         // Consistent Efield (zero current) computation.
         corrEfield_gf.ProjectCoefficient(*corrEfield_pcf);

         // Save the distribution function at a given point.
         if (right_proc_point)
		 {
            ostringstream fe_file_name;
            fe_file_name << basename << "_" << ti
                      << "_fe_point.txt";
            ofstream fe_ofs(fe_file_name.str().c_str());
            fe_ofs.precision(8);

            fe_ofs << "#x  v  f0  f1x  f0*v^2  0.5*me*f1x*v^5 0.5*me*f0*v^5\n";
            for (int i = 0; i < v_point.size(); i++)
            {
               fe_ofs << xpoint << " " << v_point[i] << " " 
                      << f0_v_point[i] << " "
                      << f1x_v_point[i] << " " << f0v2_v_point[i] << " "
                      << mehalff1xv5_v_point[i] << " "
                      << mehalff0v5_v_point[i] << endl;
            }
            fe_ofs.close();
         }
         // Save the distribution function at the max heat flux point.
         if (right_proc_pointmax)
		 {
            ostringstream fe_file_name;
            fe_file_name << basename << "_" << ti
                      << "_fe_pointmax.txt";
            ofstream fe_ofs(fe_file_name.str().c_str());
            fe_ofs.precision(8);

            fe_ofs << "#x  v  f0  f1x  f0*v^2  0.5*me*f1x*v^5 0.5*me*f0*v^5\n";
            for (int i = 0; i < v_pointmax.size(); i++)
            {
               fe_ofs << xpointmax << " " << v_pointmax[i] << " " 
                      << f0_v_pointmax[i] << " "
                      << f1x_v_pointmax[i] << " " << f0v2_v_pointmax[i] << " "
                      << mehalff1xv5_v_pointmax[i] << " "
                      << mehalff0v5_v_pointmax[i] << endl;
            }
            fe_ofs.close();
         }

         // Save spatial profiles of plasma and transport quantities.
         int NpointsPerElement = 20;
         int Nelements = x_gf.FESpace()->GetNE();
         int Npoints = NpointsPerElement * Nelements;
         double x[Npoints], rho[Npoints], Te[Npoints], j[Npoints], Ex[Npoints], 
                q[Npoints], corrE[Npoints], ne[Npoints], zbar[Npoints];
         IntegrationPoint ip;
         double dip = 1. / NpointsPerElement;
         int p = 0;
         for (int elNo = 0; elNo < Nelements; elNo++)
         {
            for (int elpoint = 0; elpoint < NpointsPerElement; elpoint++)
            {
               ip.Set3((elpoint + 0.5)*dip, 0.0, 0.0);
               x[p] = x_gf.GetValue(elNo, ip);
               rho[p] = rho_gf.GetValue(elNo, ip);
               Te[p] = e_gf.GetValue(elNo, ip);
               j[p] = jC_gf.GetValue(elNo, ip);
               Ex[p] = Efield_gf.GetValue(elNo, ip, 1);
			   q[p] = hflux_gf.GetValue(elNo, ip);
               corrE[p] = corrEfield_gf.GetValue(elNo, ip);
			   ne[p] = inElectronDens.GetValue(x[p]); 
			   zbar[p] = inZbar.GetValue(x[p]);
			   p++;
            }
         }
         ostringstream profiles_file_name;
         profiles_file_name << basename << "_" << ti << "_profiles." 
                            << setfill('0') << setw(6) << myid;
         ofstream profiles_ofs(profiles_file_name.str().c_str());
         profiles_ofs.precision(8);
         profiles_ofs << "# x  rho Te j Ex q corrE ne Zbar\n";
         for (int i = 0; i < Npoints; i++)
         {
            profiles_ofs << x[i] << " " << rho[i] << " " << Te[i] << " "
                         << j[i] << " " << Ex[i] << " " << q[i] << " "
                         << corrE[i] << " " << ne[i] << " " << zbar[i] << endl;
         }
         profiles_ofs.close();
///////////////////////////////////////////////////////////////
///// C7 nonlocal solver //////////////////////////////////////
///////////////////////////////////////////////////////////////

         if (visualization || visit || gfprint) { oper.ComputeDensity(rho_gf); }
         if (visualization)
         {
            int Wx = 0, Wy = 0; // window position
            int Ww = 350, Wh = 350; // window size
            int offx = Ww+10; // window offsets

            VisualizeField(vis_rho, vishost, visport, rho_gf,
                           "Density", Wx, Wy, Ww, Wh);
            //Wx += offx;
            //VisualizeField(vis_v, vishost, visport,
            //               v_gf, "Velocity", Wx, Wy, Ww, Wh);
            Wx += offx;
            VisualizeField(vis_j, vishost, visport, jC_gf,
                           "Current", Wx, Wy, Ww, Wh);		
            Wx += offx;
            VisualizeField(vis_e, vishost, visport, e_gf,
                           "T", Wx, Wy, Ww,Wh);

            Wx = 0;
            Wy +=offx;
            VisualizeField(vis_Efield, vishost, visport, Efield_gf, "Efield", 
                           Wx, Wy, Ww, Wh);
            //VisualizeField(vis_Kn, vishost, visport, Kn_gf, "Kn", 
            //               Wx, Wy, Ww, Wh);
            //Wx += offx;
            //VisualizeField(vis_f0, vishost, visport, intf0_gf,
            //               "int(f0 4pi v^2)dv", Wx, Wy, Ww, Wh);
            Wx += offx;
            VisualizeField(vis_hflux, vishost, visport, hflux_gf,
                           "Heat flux", Wx, Wy, Ww, Wh);
         }

         if (visit)
         {
            visit_dc.SetCycle(ti);
            visit_dc.SetTime(t);
            visit_dc.Save();
         }

         if (gfprint)
         {
            ostringstream mesh_name, rho_name, v_name, e_name;
            mesh_name << basename << "_" << ti
                      << "_mesh." << setfill('0') << setw(6) << myid;
            rho_name  << basename << "_" << ti
                      << "_rho." << setfill('0') << setw(6) << myid;
            v_name << basename << "_" << ti
                   << "_v." << setfill('0') << setw(6) << myid;
            e_name << basename << "_" << ti
                   << "_e." << setfill('0') << setw(6) << myid;

            ofstream mesh_ofs(mesh_name.str().c_str());
            mesh_ofs.precision(8);
            pmesh->Print(mesh_ofs);
            mesh_ofs.close();

            ofstream rho_ofs(rho_name.str().c_str());
            rho_ofs.precision(8);
            rho_gf.Save(rho_ofs);
            rho_ofs.close();

            ofstream v_ofs(v_name.str().c_str());
            v_ofs.precision(8);
            v_gf.Save(v_ofs);
            v_ofs.close();

            ofstream e_ofs(e_name.str().c_str());
            e_ofs.precision(8);
            e_gf.Save(e_ofs);
            e_ofs.close();

            ostringstream j_name, Kn_name, hflux_name;
            j_name << basename << "_" << ti
                   << "_j." << setfill('0') << setw(6) << myid;
            Kn_name << basename << "_" << ti
                   << "_Kn." << setfill('0') << setw(6) << myid;
            hflux_name << basename << "_" << ti
                   << "_hflux." << setfill('0') << setw(6) << myid;

            ofstream j_ofs(j_name.str().c_str());
            j_ofs.precision(8);
            jC_gf.Save(j_ofs);
            j_ofs.close();

            //ofstream Kn_ofs(Kn_name.str().c_str());
            //Kn_ofs.precision(8);
            //Kn_gf.Save(Kn_ofs);
            //Kn_ofs.close();

            ofstream hflux_ofs(hflux_name.str().c_str());
            hflux_ofs.precision(8);
            hflux_gf.Save(hflux_ofs);
            hflux_ofs.close();
         }
      }
   }

   switch (ode_solver_type)
   {
      case 2: steps *= 2; break;
      case 3: steps *= 3; break;
      case 4: steps *= 4; break;
      case 6: steps *= 6;
   }
   if (mpi.Root()) { cout << "Hydrodynamics kernel timer:" << endl << flush; }
   oper.PrintTimingData(mpi.Root(), steps);
   if (mpi.Root()) { cout << "C7 kernel timer:" << endl << flush; }
   c7oper.PrintTimingData(mpi.Root(), steps);

   if (visualization)
   {
      vis_v.close();
      vis_e.close();
   }

   // Free the used memory.
   delete ode_solver;
   delete pmesh;
   delete tensors1D; 
   delete c7ode_solver;

   return 0;
}

void display_banner(ostream & os)
{
   os << endl
      << "       __                __                 " << endl
      << "      / /   ____  ____  / /_  ____  _____   " << endl
      << "     / /   / __ `/ __ `/ __ \\/ __ \\/ ___/ " << endl
      << "    / /___/ /_/ / /_/ / / / / /_/ (__  )    " << endl
      << "   /_____/\\__,_/\\__, /_/ /_/\\____/____/  " << endl
      << "               /____/                       " << endl << endl;
}

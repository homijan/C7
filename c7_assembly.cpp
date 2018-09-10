// Copyright (c) 2017, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-734707. All Rights
// reserved. See files LICENSE and NOTICE for details.
//
// This file is part of CEED, a collection of benchmarks, miniapps, software
// libraries and APIs for efficient high-order finite element and spectral
// element discretizations for exascale applications. For more information and
// source code availability see http://github.com/ceed.
//
// The CEED research is supported by the Exascale Computing Project (17-SC-20-SC)
// a collaborative effort of two U.S. Department of Energy organizations (Office
// of Science and the National Nuclear Security Administration) responsible for
// the planning and preparation of a capable exascale ecosystem, including
// software, applications, hardware, advanced system engineering and early
// testbed platforms, in support of the nation's exascale computing imperative.

#include "c7_assembly.hpp"

#ifdef MFEM_USE_MPI

using namespace std;

namespace mfem
{

namespace nth
{

void DensityIntegrator::AssembleRHSElementVect(const FiniteElement &fe,
                                               ElementTransformation &Tr,
                                               Vector &elvect)
{
   const int ip_cnt = IntRule->GetNPoints();
   Vector shape(fe.GetDof());

   elvect.SetSize(fe.GetDof());
   elvect = 0.0;

   for (int q = 0; q < ip_cnt; q++)
   {
      fe.CalcShape(IntRule->IntPoint(q), shape);
      // Note that rhoDetJ = rho0DetJ0.
      shape *= quad_data.rho0DetJ0w(Tr.ElementNo*ip_cnt + q);
      elvect += shape;
   }
}

void Mass0Integrator::AssembleElementMatrix2(const FiniteElement &trial_fe,
                                             const FiniteElement &test_fe,
                                             ElementTransformation &Trans,
                                             DenseMatrix &elmat)
{
   const int ip_cnt = IntRule->GetNPoints();
   const int te_nd = test_fe.GetDof();
   const int tr_nd = trial_fe.GetDof();
   Vector tr_shape(tr_nd), te_shape(te_nd);

   elmat.SetSize(te_nd, tr_nd);
   elmat = 0.0;

   for (int q = 0; q < ip_cnt; q++)
   {
      test_fe.CalcShape(IntRule->IntPoint(q), te_shape);
      trial_fe.CalcShape(IntRule->IntPoint(q), tr_shape);
      te_shape *= GetIntegrator(Trans.ElementNo*ip_cnt + q);
      AddMultVWt(te_shape, tr_shape, elmat);
   }  
}

void Mass1Integrator::AssembleElementMatrix2(const FiniteElement &trial_fe,
                                             const FiniteElement &test_fe,
                                             ElementTransformation &Trans,
                                             DenseMatrix &elmat)
{
   const int nqp = IntRule->GetNPoints();
   const int dim = trial_fe.GetDim();
   const int zone_id = Trans.ElementNo;
   const int h1dofs_cnt = test_fe.GetDof();
   const int te_nd = h1dofs_cnt;
   const int tr_nd = h1dofs_cnt;

   elmat.SetSize(te_nd*dim, tr_nd*dim);
   elmat = 0.0;

   DenseMatrix partelmat(te_nd, tr_nd);
   // It is a fake AssembleMatrix2 function, since only test_fe basis is used.
   Vector shape1(te_nd);

   for (int q = 0; q < nqp; q++)
   {
      const IntegrationPoint &ip = IntRule->IntPoint(q);

      // Form stress:grad_shape at the current point.
      test_fe.CalcShape(ip, shape1);
      MultVWt(shape1, shape1, partelmat);
      partelmat *= GetIntegrator(zone_id*nqp + q);
      for (int vd = 0; vd < dim; vd++) // f1 components.
      {
         elmat.AddMatrix(partelmat, te_nd*vd, tr_nd*vd);
      }
   }
}

void DivIntegrator::AssembleElementMatrix2(const FiniteElement &trial_fe,
                                             const FiniteElement &test_fe,
                                             ElementTransformation &Trans,
                                             DenseMatrix &elmat)
{
   const int nqp = IntRule->GetNPoints();
   const int dim = trial_fe.GetDim();
   const int zone_id = Trans.ElementNo;
   const int h1dofs_cnt = test_fe.GetDof();
   const int l2dofs_cnt = trial_fe.GetDof();

   elmat.SetSize(h1dofs_cnt*dim, l2dofs_cnt);
   elmat = 0.0;

   DenseMatrix vshape(h1dofs_cnt, dim), loc_force(h1dofs_cnt, dim);
   Vector shape(l2dofs_cnt), Vloc_force(loc_force.Data(), h1dofs_cnt*dim);

   for (int q = 0; q < nqp; q++)
   {
      const IntegrationPoint &ip = IntRule->IntPoint(q);

      // Form stress:grad_shape at the current point.
      test_fe.CalcDShape(ip, vshape);
      for (int i = 0; i < h1dofs_cnt; i++)
      {
         for (int vd = 0; vd < dim; vd++) // Velocity components.
         {
            loc_force(i, vd) = 0.0;
            for (int gd = 0; gd < dim; gd++) // Gradient components.
            {
               loc_force(i, vd) += 
                  GetIntegrator(zone_id*nqp + q, vd, gd) * vshape(i,gd);
            }
         }
      }

      trial_fe.CalcShape(ip, shape);
      AddMultVWt(Vloc_force, shape, elmat);
   }
}

void VdotIntegrator::AssembleElementMatrix2(const FiniteElement &trial_fe,
                                            const FiniteElement &test_fe,
                                            ElementTransformation &Trans,
                                            DenseMatrix &elmat)
{
   const int nqp = IntRule->GetNPoints();
   const int dim = trial_fe.GetDim();
   const int zone_id = Trans.ElementNo;
   const int te_nd = test_fe.GetDof(); //h1dofs_cnt
   const int tr_nd = trial_fe.GetDof(); //l2dofs_cnt

   elmat.SetSize(te_nd*dim, tr_nd);
   elmat = 0.0;

   DenseMatrix partelmat(te_nd, tr_nd);
   Vector shape0(tr_nd), shape1(te_nd);

   for (int q = 0; q < nqp; q++)
   {
      const IntegrationPoint &ip = IntRule->IntPoint(q);

      // Form stress:grad_shape at the current point.
      test_fe.CalcShape(ip, shape1);
      trial_fe.CalcShape(ip, shape0);
      MultVWt(shape1, shape0, partelmat);

      for (int vd = 0; vd < dim; vd++) // f1 components.
      {
         elmat.AddMatrix(GetIntegrator(zone_id*nqp + q, vd), partelmat,
                         te_nd*vd, 0);
      }
   }
}

void VcrossIntegrator::AssembleElementMatrix2(const FiniteElement &trial_fe,
                                              const FiniteElement &test_fe,
                                              ElementTransformation &Trans,
                                              DenseMatrix &elmat)
{
   const int nqp = IntRule->GetNPoints();
   const int dim = trial_fe.GetDim();
   const int zone_id = Trans.ElementNo;
   const int h1dofs_cnt = test_fe.GetDof();

   elmat.SetSize(h1dofs_cnt*dim, h1dofs_cnt*dim);
   elmat = 0.0;

   DenseMatrix loc_w(h1dofs_cnt, dim), loc_vxw(h1dofs_cnt, dim);
   // It is a fake AssembleMatrix2 function, since only test_fe basis is used.
   Vector shape1(h1dofs_cnt);
   Vector Vloc_w(loc_w.Data(), h1dofs_cnt*dim),
          Vloc_vxw(loc_vxw.Data(), h1dofs_cnt*dim);

   for (int q = 0; q < nqp; q++)
   {
      const IntegrationPoint &ip = IntRule->IntPoint(q);

      // Form stress:grad_shape at the current point.
      test_fe.CalcShape(ip, shape1);
      for (int i = 0; i < h1dofs_cnt; i++)
      {
         for (int vd = 0; vd < dim; vd++) // f1 components.
         {
            loc_w(i, vd) = shape1(i);
            loc_vxw(i, vd) = GetIntegrator(zone_id*nqp + q, vd) * shape1(i);
         }
      }

      AddMultVWt(Vloc_w, Vloc_vxw, elmat);
   }
}

double Mass0NuIntegrator::GetIntegrator(int i)
{
   return quad_data.nu_invrho(i) * quad_data.rho0DetJ0w(i);
}

double invMass0NuIntegrator::GetIntegrator(int i)
{
   return quad_data.nu_invrho(i) * quad_data.rho0DetJ0w(i);
}

double Mass1NuIntegrator::GetIntegrator(int i)
{
   return quad_data.nu_invrho(i) * quad_data.rho0DetJ0w(i);
}

double Mass1NutIntegrator::GetIntegrator(int i)
{
   return quad_data.nut_invrho(i) * quad_data.rho0DetJ0w(i);
}

double Divf0Integrator::GetIntegrator(int i, int vd, int gd)
{
   return quad_data.stress0JinvT(vd)(i, gd);
}

double Divf1Integrator::GetIntegrator(int i, int vd, int gd)
{
   return quad_data.stress1JinvT(vd)(i, gd);
}

double EfieldScatterIntegrator::GetIntegrator(int i, int vd)
{
   return quad_data.Escatter_invrho(i, vd) * quad_data.rho0DetJ0w(i);
}

double EfieldScaledIntegrator::GetIntegrator(int i, int vd)
{
   return quad_data.Escaled_invrho(i, vd) * quad_data.rho0DetJ0w(i);
}

double AEfieldScaledIntegrator::GetIntegrator(int i, int vd)
{
   return quad_data.AEscaled_invrho(i, vd) * quad_data.rho0DetJ0w(i);
}

double AEfieldIntegrator::GetIntegrator(int i, int vd)
{
   return quad_data.AE_invrho(i, vd) * quad_data.rho0DetJ0w(i);
}

double BfieldIntegrator::GetIntegrator(int i, int vd)
{
   return quad_data.B_invrho(i, vd) * quad_data.rho0DetJ0w(i);
}

} // namespace nth

} // namespace mfem

#endif // MFEM_USE_MPI

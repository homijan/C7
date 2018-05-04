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

#ifndef MFEM_NTH_C7_ASSEMBLY
#define MFEM_NTH_C7_ASSEMBLY

#include "mfem.hpp"


#ifdef MFEM_USE_MPI

#include <memory>
#include <iostream>

namespace mfem
{

namespace nth
{

// Container for all data needed at quadrature points.
struct QuadratureData
{
   // TODO: use QuadratureFunctions?

   // Reference to physical Jacobian for the initial mesh. These are computed
   // only at time zero and stored here.
   DenseTensor Jac0inv;

   // Quadrature data used for full/partial assembly of the force operator. At
   // each quadrature point, it combines the stress, inverse Jacobian,
   // determinant of the Jacobian and the integration weight. It must be
   // recomputed in every time step.
   DenseTensor stress1JinvT, stress0JinvT;

   // Quadrature data used for full/partial assembly of the mass matrices. At
   // time zero, we compute and store (rho0 * det(J0) * qp_weight) at each
   // quadrature point. Note the at any other time, we can compute
   // rho = rho0 * det(J0) / det(J), representing the notion of pointwise mass
   // conservation.
   Vector rho0DetJ0w;

   // Mass integrators.
   Vector nuinvrho, nutinvrho, nuEinvrho;

   // The pointwise equality rho * detJ = rho0 * detJ0 is used by integrators.
   // Electric and magnetic fields. 
   DenseMatrix Escaled_invrho, AEscaled_invrho, AEinvrho, AIEinvrho, Binvrho;
   // Explicit zero moment "mass" integrator.
   Vector Ef1invvf0rho;

   // Initial length scale. This represents a notion of local mesh size. We
   // assume that all initial zones have similar size.
   double h0;

   // Estimate of the minimum time step over all quadrature points. This is
   // recomputed at every time step to achieve adaptive time stepping.
   double dt_est;

   QuadratureData(int dim, int nzones, int quads_per_zone)
      : Jac0inv(dim, dim, nzones * quads_per_zone),
        stress1JinvT(nzones * quads_per_zone, dim, dim),
        stress0JinvT(nzones * quads_per_zone, dim, dim),
        rho0DetJ0w(nzones * quads_per_zone),
        nuinvrho(nzones * quads_per_zone),
        nutinvrho(nzones * quads_per_zone),
        nuEinvrho(nzones * quads_per_zone),
        Escaled_invrho(nzones * quads_per_zone, dim),
        AEscaled_invrho(nzones * quads_per_zone, dim),
        AEinvrho(nzones * quads_per_zone, dim),
        AIEinvrho(nzones * quads_per_zone, dim),
        Binvrho(nzones * quads_per_zone, dim),
        Ef1invvf0rho(nzones * quads_per_zone) { }
};

// This class is used only for visualization. It assembles (rho, phi) in each
// zone, which is used by LagrangianHydroOperator::ComputeDensity to do an L2
// projection of the density.
class DensityIntegrator : public LinearFormIntegrator
{
private:
   const QuadratureData &quad_data;

public:
   DensityIntegrator(QuadratureData &quad_data_) : quad_data(quad_data_) { }

   virtual void AssembleRHSElementVect(const FiniteElement &fe,
                                       ElementTransformation &Tr,
                                       Vector &elvect);
};

// Assembles element contributions to the global velocity force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class Mass0Integrator : public BilinearFormIntegrator
{
protected:
   const QuadratureData &quad_data;

public:
   Mass0Integrator(QuadratureData &quad_data_) : quad_data(quad_data_) { }

   virtual void AssembleElementMatrix(const FiniteElement &fe,
                                      ElementTransformation &Trans,
                                      DenseMatrix &elmat)
   { AssembleElementMatrix2(fe, fe, Trans, elmat); }
   virtual void AssembleElementMatrix2(const FiniteElement &trial_fe,
                                       const FiniteElement &test_fe,
                                       ElementTransformation &Trans,
                                       DenseMatrix &elmat);
   virtual double GetIntegrator(int i) = 0;
};

// Assembles element contributions to the global velocity force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class invMass0Integrator : public Mass0Integrator
{
public:
   invMass0Integrator(QuadratureData &quad_data_) 
      : Mass0Integrator(quad_data_)  { }

   virtual void AssembleElementMatrix2(const FiniteElement &trial_fe,
                                       const FiniteElement &test_fe,
                                       ElementTransformation &Trans,
                                       DenseMatrix &elmat)
   {  
      // Compute the element Mass0 matrix.
      Mass0Integrator::AssembleElementMatrix2(trial_fe, test_fe, Trans, elmat);
      // Invert the element matrix.
      elmat.Invert();   
   } 
};

// Assembles element contributions to the global velocity force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class Mass1Integrator : public BilinearFormIntegrator
{
protected:
   const QuadratureData &quad_data;

public:
   Mass1Integrator(QuadratureData &quad_data_) : quad_data(quad_data_) { }

   virtual void AssembleElementMatrix(const FiniteElement &fe,
                                      ElementTransformation &Trans,
                                      DenseMatrix &elmat)
   { AssembleElementMatrix2(fe, fe, Trans, elmat); }
   virtual void AssembleElementMatrix2(const FiniteElement &trial_fe,
                                       const FiniteElement &test_fe,
                                       ElementTransformation &Trans,
                                       DenseMatrix &elmat);
   virtual double GetIntegrator(int i) = 0;
};

// Assembles element contributions to the global velocity force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class DivIntegrator : public BilinearFormIntegrator
{
protected:
   const QuadratureData &quad_data;

public:
   DivIntegrator(QuadratureData &quad_data_) : quad_data(quad_data_) { }

   virtual void AssembleElementMatrix2(const FiniteElement &trial_fe,
                                       const FiniteElement &test_fe,
                                       ElementTransformation &Trans,
                                       DenseMatrix &elmat);
   virtual double GetIntegrator(int q, int vd, int gd) = 0;
};

// Assembles element contributions to the global velocity force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class VdotIntegrator : public BilinearFormIntegrator
{
protected:
   const QuadratureData &quad_data;

public:
   VdotIntegrator(QuadratureData &quad_data_) : quad_data(quad_data_) { }

   virtual void AssembleElementMatrix2(const FiniteElement &trial_fe,
                                       const FiniteElement &test_fe,
                                       ElementTransformation &Trans,
                                       DenseMatrix &elmat);
   virtual double GetIntegrator(int q, int vd) = 0;
};

// Assembles element contributions to the global velocity force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class VcrossIntegrator : public BilinearFormIntegrator
{
protected:
   const QuadratureData &quad_data;

public:
   VcrossIntegrator(QuadratureData &quad_data_) : quad_data(quad_data_) { }

   virtual void AssembleElementMatrix(const FiniteElement &fe,
                                      ElementTransformation &Trans,
                                      DenseMatrix &elmat)
   { AssembleElementMatrix2(fe, fe, Trans, elmat); }
   virtual void AssembleElementMatrix2(const FiniteElement &trial_fe,
                                       const FiniteElement &test_fe,
                                       ElementTransformation &Trans,
                                       DenseMatrix &elmat);
   virtual double GetIntegrator(int q, int vd) = 0;
};

// Assembles element contributions to the global velocity force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class Mass0NuIntegrator : public Mass0Integrator
{
private:
public:
   Mass0NuIntegrator(QuadratureData &quad_data_) :
      Mass0Integrator(quad_data_) { }

   double GetIntegrator(int i);
};

// Assembles element contributions to the global velocity force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class ExplMass0Integrator : public Mass0Integrator
{
private:
public:
   ExplMass0Integrator(QuadratureData &quad_data_) :
      Mass0Integrator(quad_data_) { }

   double GetIntegrator(int i);
};

// Assembles element contributions to the global velocity force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class invMass0NuIntegrator : public invMass0Integrator
{
private:
public:
   invMass0NuIntegrator(QuadratureData &quad_data_) :
      invMass0Integrator(quad_data_) { }

   double GetIntegrator(int i);
};

// Assembles element contributions to the global velocity force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class invMass0NuEIntegrator : public invMass0Integrator
{
private:
public:
   invMass0NuEIntegrator(QuadratureData &quad_data_) :
      invMass0Integrator(quad_data_) { }

   double GetIntegrator(int i);
};

// Assembles element contributions to the global velocity force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class Mass1NuIntegrator : public Mass1Integrator
{
private:
public:
   Mass1NuIntegrator(QuadratureData &quad_data_) :
      Mass1Integrator(quad_data_) { }

   double GetIntegrator(int i);
};

// Assembles element contributions to the global velocity force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class Mass1NutIntegrator : public Mass1Integrator
{
private:
public:
   Mass1NutIntegrator(QuadratureData &quad_data_) :
      Mass1Integrator(quad_data_) { }

   double GetIntegrator(int i);
};

// Assembles element contributions to the global temperature force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class Divf0Integrator : public DivIntegrator
{
private:

public:
   Divf0Integrator(QuadratureData &quad_data_) : DivIntegrator(quad_data_) { }

   double GetIntegrator(int q, int vd, int gd);
};

// Assembles element contributions to the global velocity force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class Divf1Integrator : public DivIntegrator
{
private:
public:
   Divf1Integrator(QuadratureData &quad_data_) : DivIntegrator(quad_data_) { }

   double GetIntegrator(int q, int vd, int gd);
};

// Assembles element contributions to the global temperature force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class EfieldScIntegrator : public VdotIntegrator
{
private:

public:
   EfieldScIntegrator(QuadratureData &quad_data_) : VdotIntegrator(quad_data_) 
   { }

   double GetIntegrator(int q, int vd);
};

// Assembles element contributions to the global temperature force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class AEfieldScIntegrator : public VdotIntegrator
{
private:

public:
   AEfieldScIntegrator(QuadratureData &quad_data_) :
      VdotIntegrator(quad_data_) { }

   double GetIntegrator(int q, int vd);
};

// Assembles element contributions to the global temperature force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class AEfieldIntegrator : public VdotIntegrator
{
private:

public:
   AEfieldIntegrator(QuadratureData &quad_data_) :
      VdotIntegrator(quad_data_) { }

   double GetIntegrator(int q, int vd);
};

// Assembles element contributions to the global temperature force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class AIEfieldIntegrator : public VdotIntegrator
{
private:

public:
   AIEfieldIntegrator(QuadratureData &quad_data_) :
      VdotIntegrator(quad_data_) { }

   double GetIntegrator(int q, int vd);
};

// Assembles element contributions to the global temperature force matrix.
// This class is used for the full assembly case; it's not used with partial
// assembly.
class BfieldIntegrator : public VcrossIntegrator
{
private:

public:
   BfieldIntegrator(QuadratureData &quad_data_) :
      VcrossIntegrator(quad_data_) { }

   double GetIntegrator(int q, int vd);
};

} // namespace nth

} // namespace mfem

#endif // MFEM_USE_MPI

#endif // MFEM_NTH_C7_ASSEMBLY

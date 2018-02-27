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

#ifndef MFEM_NTH_AWBSPHYSICS
#define MFEM_NTH_AWBSPHYSICS

#include "mfem.hpp"
#include "eos.hpp"

#ifdef MFEM_USE_MPI

//#include <memory>
//#include <iostream>
//#include <fstream>

namespace mfem
{

namespace nth
{

// NTH hydro coefficient including velocity dependence.
class NTHvHydroCoefficient : public HydroCoefficient
{
   void SetVelocityScale(double N_x_, double Tmax)
      { N_x_vTmax = N_x * eos->GetvTe(Tmax); }
protected:
   // Velocity is always scaled wit respect to maximum thermal velocity
   // (its multiple) so it is in (0, 1)
   double N_x, Tmax, N_x_vTmax;
   // Current particle velocity from the velocity spectra.
   double velocity;
public:
   NTHvHydroCoefficient(ParGridFunction &rho_, ParGridFunction &T_,
                        ParGridFunction &v_, Coefficient *material_, EOS *eos_)
      : HydroCoefficient(rho_, T_, v_, material_, eos_)
	  { N_x = 1.0; Tmax = 1.0; SetVelocityScale(N_x, Tmax); }
   virtual double Eval(ElementTransformation &T,
      const IntegrationPoint &ip)
      { double rho = rho_gf.GetValue(T.ElementNo, ip); Eval(T, ip, rho); }
   virtual double Eval(ElementTransformation &T,
      const IntegrationPoint &ip, double rho) = 0; 
   void SetVelocity(double v_) { velocity = v_; }
   void SetThermalVelocityMultiple(double N_x_)
      { N_x = N_x_; SetVelocityScale(N_x, Tmax); }
   void SetTmax(double Tmax_)
      { Tmax = Tmax_; SetVelocityScale(N_x, Tmax); }
   double GetVelocityScale() { return N_x_vTmax; }
   double GetRho(ElementTransformation &T, const IntegrationPoint &ip)
      { return rho_gf.GetValue(T.ElementNo, ip); }
};

// Master of physics for nonlocal transport.
class AWBSMasterOfPhysics
{
protected:
public:
   // members
   NTHvHydroCoefficient *mspei_pcf, *mspee_pcf, *sourceF0_pcf;
   VectorCoefficient *Efield_pcf, *Bfield_pcf;
   // methods
   AWBSMasterOfPhysics(NTHvHydroCoefficient *mspei_, 
                       NTHvHydroCoefficient *mspee_,
                       NTHvHydroCoefficient *source_, 
                       VectorCoefficient *Efield_,
                       VectorCoefficient *Bfield_)
      : mspei_pcf(mspei_),  mspee_pcf(mspee_), sourceF0_pcf(source_),
        Efield_pcf(Efield_), Bfield_pcf(Bfield_) { }

   void SetThermalVelocityMultiple(double vTmultiple)
   { 
      mspei_pcf->SetThermalVelocityMultiple(vTmultiple);
      mspee_pcf->SetThermalVelocityMultiple(vTmultiple);
      sourceF0_pcf->SetThermalVelocityMultiple(vTmultiple);
   }
   void SetTmax(double glob_Tmax)
   { 
      mspei_pcf->SetTmax(glob_Tmax); 
      mspee_pcf->SetTmax(glob_Tmax);
      sourceF0_pcf->SetTmax(glob_Tmax);
   }

   ~AWBSMasterOfPhysics() { }
};

// Classical mean-stopping-power (electron-ion) coefficient.
class ClassicalMeanStoppingPower : public NTHvHydroCoefficient
{
protected:
public:
   ClassicalMeanStoppingPower(ParGridFunction &rho_, ParGridFunction &Te_,
                              ParGridFunction &v_, Coefficient *material_,
                              EOS *eos_)
      : NTHvHydroCoefficient(rho_, Te_, v_, material_, eos_) {}
   //virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip);
   virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip,
                       double rho);
};

// Classical electron-ion mean-stopping-power coefficient.
class ClassicalAWBSMeanStoppingPower : public ClassicalMeanStoppingPower
{
protected:
public:
   ClassicalAWBSMeanStoppingPower(ParGridFunction &rho_, 
                                  ParGridFunction &Te_,
                                  ParGridFunction &v_, 
                                  Coefficient *material_,
                                  EOS *eos_)
      : ClassicalMeanStoppingPower(rho_, Te_, v_, material_, eos_) {}
   virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip,
                       double rho);
};

// General mean-free-path coefficient.
class MeanFreePath
{
protected:
public:
   MeanFreePath() {}
   virtual double EvalThermalMFP(ElementTransformation &T, 
                                 const IntegrationPoint &ip, double rho) = 0;
};

// Classical mean-free-path coefficient.
class ClassicalMeanFreePath : public MeanFreePath, 
                              public ClassicalMeanStoppingPower
{
protected:
public:
   ClassicalMeanFreePath(ParGridFunction &rho_, ParGridFunction &Te_,
                         ParGridFunction &v_, Coefficient *material_, 
                         EOS *eos_)
      : ClassicalMeanStoppingPower(rho_, Te_, v_, material_, eos_) {}
   virtual double EvalThermalMFP(ElementTransformation &T, 
                                 const IntegrationPoint &ip, double rho);
};

// Classical Kn(mean-stopping-power) coefficient.
class KnudsenNumber : public HydroCoefficient
{
protected:
   MeanFreePath *mfp;
public:
   KnudsenNumber(ParGridFunction &rho_, ParGridFunction &Te_,
                 ParGridFunction &v_, Coefficient *material_,
                 EOS *eos_, MeanFreePath *mfp_)
      : HydroCoefficient(rho_, Te_, v_, material_, eos_), mfp(mfp_) {}
   virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip);
};

// Classical Lorentz approximation E field coefficient.
class LorentzEfield : public HydroCoefficient, public VectorCoefficient
{
protected:
   double S0; 
public:
   LorentzEfield(int dim_, ParGridFunction &rho_, ParGridFunction &Te_,
                 ParGridFunction &v_, Coefficient *material_, EOS *eos_)
      : HydroCoefficient(rho_, Te_, v_, material_, eos_), 
        VectorCoefficient(dim_) { S0 = 1.0; }
   virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip);
   virtual void Eval(Vector &V, ElementTransformation &T,
                     const IntegrationPoint &ip);
   void SetScale0(double S0_) { S0 = S0_; }
};

// AWBS source coefficient.
class AWBSF0Source : public NTHvHydroCoefficient
{
protected:
   double S0;
public:
   AWBSF0Source(ParGridFunction &rho_, ParGridFunction &Te_,
                ParGridFunction &v_, Coefficient *material_, EOS *eos_)
      : NTHvHydroCoefficient(rho_, Te_, v_, material_, eos_) { S0 = 1.0; }
   double Eval(ElementTransformation &T, const IntegrationPoint &ip);
   double Eval(ElementTransformation &T, const IntegrationPoint &ip,
               double rho);
   void SetScale0(double S0_) { S0 = S0_; }
};

extern double sigma;
extern double coulLog; // TMP will be moved to the eos.

} // namespace nth

} // namespace mfem

#endif // MFEM_USE_MPI

#endif // MFEM_NTH_AWBSPHYSICS

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
   //void SetVelocityScale(double N_x_, double Tmax)
   //   { N_x_vTmax = N_x * eos->GetvTe(Tmax); }
protected:
   // Velocity is always scaled wit respect to maximum thermal velocity
   // (its multiple) so it is in (0, 1)
   //double N_x, Tmax, N_x_vTmax;
   // Current particle velocity from the velocity spectra.
   double velocity_real;
public:
   NTHvHydroCoefficient(ParGridFunction &rho_, ParGridFunction &T_,
                        ParGridFunction &v_, Coefficient *material_, EOS *eos_)
      : HydroCoefficient(rho_, T_, v_, material_, eos_) { }
	  //{ N_x = 1.0; Tmax = 1.0; SetVelocityScale(N_x, Tmax); }
   virtual double Eval(ElementTransformation &T,
      const IntegrationPoint &ip)
      { double rho = rho_gf.GetValue(T.ElementNo, ip); Eval(T, ip, rho); }
   virtual double Eval(ElementTransformation &T,
      const IntegrationPoint &ip, double rho) = 0; 
   void SetVelocityReal(double v_) { velocity_real = v_; }
   //void SetThermalVelocityMultiple(double N_x_)
   //   { N_x = N_x_; SetVelocityScale(N_x, Tmax); }
   //void SetTmax(double Tmax_)
   //   { Tmax = Tmax_; SetVelocityScale(N_x, Tmax); }
   //double GetVelocityMultiple() { return N_x; }
   //double GetVelocityScale() { return N_x_vTmax; }
   double GetRho(ElementTransformation &T, const IntegrationPoint &ip)
      { return rho_gf.GetValue(T.ElementNo, ip); }
   double GetTe(ElementTransformation &T, const IntegrationPoint &ip)
      { return Te_gf.GetValue(T.ElementNo, ip); }
   double GetvTe(ElementTransformation &T, const IntegrationPoint &ip)
      { return eos->GetvTe(GetTe(T, ip)); }
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
   double corrAWBS;
public:
   ClassicalAWBSMeanStoppingPower(ParGridFunction &rho_, 
                                  ParGridFunction &Te_,
                                  ParGridFunction &v_, 
                                  Coefficient *material_,
                                  EOS *eos_)
      : ClassicalMeanStoppingPower(rho_, Te_, v_, material_, eos_) 
      { corrAWBS = 1.0; }
   virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip,
                       double rho);
   virtual void SetCorrAWBS(double corrAWBS_) { corrAWBS = corrAWBS_; }
};

/*
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
*/

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

// General entity operating on kinetic distribution function.
class GeneralKineticCoefficient : public Coefficient, public VectorCoefficient
{
protected:
   double velocity_real;
   // GridFunctions related to F0 and F1, e.g. F0, F1, dF0dv, dF1dv, etc.
   ParGridFunction *F0, *F1, *dF0, *dF1;
   NTHvHydroCoefficient *mspei_pcf, *mspee_pcf;
   VectorCoefficient *Efield_pcf;
public:
   GeneralKineticCoefficient(int dim_, NTHvHydroCoefficient *mspei_pcf_,
                             NTHvHydroCoefficient *mspee_pcf_)
      : F0(NULL), F1(NULL), VectorCoefficient(dim_), mspei_pcf(mspei_pcf_), 
      mspee_pcf(mspee_pcf_) { }

   virtual void SetF0(ParGridFunction *F0_) { F0 = F0_; }
   virtual void SetF1(ParGridFunction *F1_) { F1 = F1_; }
   virtual void SetdF0(ParGridFunction *dF0_) { dF0 = dF0_; }
   virtual void SetdF1(ParGridFunction *dF1_) { dF1 = dF1_; }
   virtual void SetEfield(VectorCoefficient *Efield_) { Efield_pcf = Efield_; }
   virtual void SetVelocityReal(double v_) { velocity_real = v_; }

   virtual double Eval(ElementTransformation &T, 
                       const IntegrationPoint &ip) = 0;
   virtual void Eval(Vector &V, ElementTransformation &T,
                     const IntegrationPoint &ip) = 0;
};

// P1 Efield Lorentz force kinetic coefficient.
class P1a0KineticCoefficient : public GeneralKineticCoefficient
{
protected:
public:
   P1a0KineticCoefficient(int dim_, NTHvHydroCoefficient *mspei_pcf_,
                          NTHvHydroCoefficient *mspee_pcf_)
      : GeneralKineticCoefficient(dim_, mspei_pcf_, mspee_pcf_) { }
   virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip);
   virtual void Eval(Vector &V, ElementTransformation &T,
                     const IntegrationPoint &ip) 
   { std::cout << "Not defined..." << std::endl; }
};

// P1 background force kinetic coefficient.
class P1b0KineticCoefficient : public GeneralKineticCoefficient
{
protected:
public:
   P1b0KineticCoefficient(int dim_, NTHvHydroCoefficient *mspei_pcf_,
                          NTHvHydroCoefficient *mspee_pcf_)
      : GeneralKineticCoefficient(dim_, mspei_pcf_, mspee_pcf_) { }
   virtual double Eval(ElementTransformation &T, 
                       const IntegrationPoint &ip)
   { std::cout << "Not defined..." << std::endl; }
   virtual void Eval(Vector &V, ElementTransformation &T,
                     const IntegrationPoint &ip);
};

// P1 magnetic force kinetic coefficient.
class P1b1KineticCoefficient : public GeneralKineticCoefficient
{
protected:
public:
   P1b1KineticCoefficient(int dim_, NTHvHydroCoefficient *mspei_pcf_,
                          NTHvHydroCoefficient *mspee_pcf_)
      : GeneralKineticCoefficient(dim_, mspei_pcf_, mspee_pcf_) { }
   virtual double Eval(ElementTransformation &T, 
                       const IntegrationPoint &ip)
   { std::cout << "Not defined..." << std::endl; }
   virtual void Eval(Vector &V, ElementTransformation &T,
                     const IntegrationPoint &ip);
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

// Master of physics for nonlocal transport.
class AWBSMasterOfPhysics
{
protected:
   // Integration over velocity is always scaled wit respect to 
   // a global maximum of thermal velocity (its multiple) so 
   // velocity spans (0, 1).
   double N_x, Tmax, N_x_vTmax;
   // General equation of state.
   EOS *eos;
public:
   // members
   NTHvHydroCoefficient *mspei_pcf, *mspee_pcf, *sourceF0_pcf;
   VectorCoefficient *Efield_pcf, *Bfield_pcf;
   P1a0KineticCoefficient *P1a0_pcf;
   P1b0KineticCoefficient *P1b0_pcf;
   P1b1KineticCoefficient *P1b1_pcf;
   // methods
   AWBSMasterOfPhysics(int dim_, NTHvHydroCoefficient *mspei_, 
                       NTHvHydroCoefficient *mspee_,
                       NTHvHydroCoefficient *source_, 
                       VectorCoefficient *Efield_,
                       VectorCoefficient *Bfield_, EOS *eos_)
      : mspei_pcf(mspei_),  mspee_pcf(mspee_), sourceF0_pcf(source_), 
        Efield_pcf(Efield_), Bfield_pcf(Bfield_), eos(eos_)
        { 
           P1a0_pcf = new P1a0KineticCoefficient(dim_, mspei_, mspee_); 
           P1b0_pcf = new P1b0KineticCoefficient(dim_, mspei_, mspee_); 
		   P1b1_pcf = new P1b1KineticCoefficient(dim_, mspei_, mspee_);
		   P1a0_pcf->SetEfield(Efield_);
        }
   void SetVelocityScale(double N_x_, double Tmax)
      { N_x = N_x_; N_x_vTmax = N_x_ * eos->GetvTe(Tmax); }

   //void SetThermalVelocityMultiple(double vTmultiple)
   //{ 
   //   N_x = vTmultiple;
   //   mspei_pcf->SetThermalVelocityMultiple(vTmultiple);
   //   mspee_pcf->SetThermalVelocityMultiple(vTmultiple);
   //   sourceF0_pcf->SetThermalVelocityMultiple(vTmultiple);
   //}
   //void SetTmax(double glob_Tmax)
   //{ 
   //   mspei_pcf->SetTmax(glob_Tmax); 
   //   mspee_pcf->SetTmax(glob_Tmax);
   //   sourceF0_pcf->SetTmax(glob_Tmax);
   //}
   double GetVelocityMultiple() { return N_x; }
   double GetVelocityScale() { return N_x_vTmax; }

   ~AWBSMasterOfPhysics() 
   { 
      delete P1a0_pcf;
      delete P1b0_pcf;
      delete P1b1_pcf;
   }
};

// Post-processing coefficient of generalized Ohm's law current.
class OhmCurrentCoefficient : public VectorCoefficient
{
protected:
   ParGridFunction *a0_pgf, *b0_pgf, *b1_pgf; 
   VectorCoefficient *Efield_pcf, *Bfield_pcf;
public:
   OhmCurrentCoefficient(int dim_, ParGridFunction *a0_pgf_, 
                         ParGridFunction *b0_pgf_, ParGridFunction *b1_pgf_,
                         VectorCoefficient *Efield_pcf_,
                         VectorCoefficient *Bfield_pcf_)
      : VectorCoefficient(dim_), a0_pgf(a0_pgf_), b0_pgf(b0_pgf_), 
      b1_pgf(b1_pgf_), Efield_pcf(Efield_pcf_), Bfield_pcf(Bfield_pcf_) { }

   virtual void Eval(Vector &V, ElementTransformation &T,
                     const IntegrationPoint &ip);
};

// Post-processing coefficient of generalized Ohm's law E field 
// corresponding to reduced current.
class OhmEfieldCoefficient : public VectorCoefficient
{
protected:
   double alpha;
   ParGridFunction *jC_pgf, *a0_pgf, *b0_pgf, *b1_pgf; 
   VectorCoefficient *Bfield_pcf;
public:
   OhmEfieldCoefficient(int dim_, ParGridFunction *jC_pgf_,
                        ParGridFunction *a0_pgf_, ParGridFunction *b0_pgf_, 
                        ParGridFunction *b1_pgf_,
                        VectorCoefficient *Bfield_pcf_)
      : VectorCoefficient(dim_), jC_pgf(jC_pgf_), a0_pgf(a0_pgf_), 
      b0_pgf(b0_pgf_), b1_pgf(b1_pgf_), Bfield_pcf(Bfield_pcf_) { alpha = 0.5; }

   virtual void Eval(Vector &V, ElementTransformation &T,
                     const IntegrationPoint &ip);
   virtual void SetAlpha(double alpha_) { alpha = alpha_; }
};

extern double sigma;
extern double coulLog; // TMP will be moved to the eos.

} // namespace nth

} // namespace mfem

#endif // MFEM_USE_MPI

#endif // MFEM_NTH_AWBSPHYSICS

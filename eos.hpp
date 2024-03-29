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

#ifndef MFEM_NTH_EOS
#define MFEM_NTH_EOS

#include "mfem.hpp"

#ifdef MFEM_USE_MPI

namespace mfem
{

namespace nth
{

// Generic hydro equation of state (EOS) class,
// providing any physics related evaluation needed in NTH.
class EOS
{
protected:
   // Fundamental constants of nature.
   double kB, c, hbar, G;
   // Corresponding masses of electron and proton.
   double me;
public:
   EOS(double me_ = 1.0, double kB_ = 1.0, double c_ = 1.0,
       double hbar_ = 1.0, double G_ = 1.0)
      { kB = kB_, c = c_, hbar = hbar_, G = G_; me = me_; }
   // Get fundamental physical constants.
   double GetkB() {return kB; }
   double Getc() {return c; }
   double Getkhbar() {return hbar; }
   double GetG() {return G; }
   double Getme() {return me; }
   // Thermal velocity.
   double GetvTe(double Te) { return sqrt(kB * Te / me); }
   // Get Thermodynamic values.
   virtual double GetIonDensity(double index, double rho) = 0;
   virtual double GetElectronDensity(double index, double rho, double Te)
           { return GetZbar(index, rho, Te) * GetIonDensity(index, rho); }
   virtual double GetZbar(double index, double rho, double Te) = 0;
   virtual double GetPe(double index, double rho, double Te) = 0;
   virtual double GetPi(double index, double rho, double Te) = 0;
   virtual double GetEe(double index, double rho, double Te) = 0;
   virtual double GetEi(double index, double rho, double Te) = 0;
   virtual double GetdEedTe(double index, double rho, double Te) = 0;
   virtual double GetdEidTi(double index, double rho, double Te) = 0;
   virtual double GetdEedrho(double index, double rho, double Te) = 0;
   virtual double GetdEidrho(double index, double rho, double Te) = 0;
   virtual double GetSoundSpeed(double index, double rho, double Te) = 0;

   ~EOS() { }
};

// The simplest equation of state to be used in hydrodynamic simulations. 
// Usually it is referred to as ideal gas.
class IGEOS : public EOS
{
protected:
   // Mean ionization, here defined as constant.
   double Zbar;
   // Ion mass [g], here defined as constant.
   double mi;
public:
   IGEOS(double me_, double kB_) : EOS(me_, kB_, 1.0, 1.0, 1.0) 
   { Zbar = 1.0; mi = 1.0; }
   // Get Thermodynamic values.
   virtual double GetIonDensity(double index, double rho) { return rho / mi; }
   virtual double GetZbar(double index, double rho, double Te) { return Zbar; }
   virtual double GetPe(double index, double rho, double Te) {}
   virtual double GetPi(double index, double rho, double Te) {}
   virtual double GetEe(double index, double rho, double Te) {}
   virtual double GetEi(double index, double rho, double Te) {}
   virtual double GetdEedTe(double index, double rho, double Te) {}
   virtual double GetdEidTi(double index, double rho, double Te) {}
   virtual double GetdEedrho(double index, double rho, double Te) {}
   virtual double GetdEidrho(double index, double rho, double Te) {}
   virtual double GetSoundSpeed(double index, double rho, double Te) {}
   // IG specific functions.
   void SetZbar(double Zbar_) { Zbar = Zbar_; }
   void SetIonMass(double mi_) { mi = mi_; }
};

// Generic hydro coefficient.
class HydroCoefficient : public Coefficient
{
protected:
   // Fluid quantities used in calculations of physics.
   ParGridFunction &rho_gf, &Te_gf, &v_gf;
   // Space dependent material coefficient.
   Coefficient *material_pcf;
   // General equation of state.
   EOS *eos;
public:
   HydroCoefficient(ParGridFunction &rho_, ParGridFunction &Te_,
                    ParGridFunction &v_, Coefficient *material_, EOS *eos_)
      : rho_gf(rho_), Te_gf(Te_), v_gf(v_), material_pcf(material_), eos(eos_)
      { }
   virtual double Eval(ElementTransformation &T,
      const IntegrationPoint &ip) = 0;
   virtual void SetEOS(EOS *eos_) { eos = eos_; }

   virtual ~HydroCoefficient() {};
};

} // namespace nth

} // namespace mfem

#endif // MFEM_USE_MPI

#endif // MFEM_NTH_EOS

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

#include "AWBSphysics.hpp"

#ifdef MFEM_USE_MPI

using namespace std;

namespace mfem
{

namespace nth
{

double sigma = 8.1027575e17; // Matching the SH diffusive flux.
double coulLog = 10.0; // TMP will be moved to the eos.

double ClassicalMeanStoppingPower::Eval(ElementTransformation &T,
                                        const IntegrationPoint &ip, double rho)
{
   // "velocity_real" must be set previously.
   double Te = Te_gf.GetValue(T.ElementNo, ip);
   double index = material_pcf->Eval(T, ip);
   double Zbar = eos->GetZbar(index, rho, Te);
   double ni = eos->GetIonDensity(index, rho, Te);
   // The ei collision frequency has standard 1/v^3 dependence, 
   // the sigma cross section given by the model sigma * ni 
   // and Zbar increases the effect of Coulomb potential in "on ion collisions".
   double nu = Zbar * Zbar * sigma * coulLog * ni / pow(velocity_real, 3.0);

   return nu;
}

double ClassicalAWBSMeanStoppingPower::Eval(ElementTransformation &T,
                                            const IntegrationPoint &ip, 
                                            double rho)
{
   double Te = Te_gf.GetValue(T.ElementNo, ip);
   double index = material_pcf->Eval(T, ip);
   double Zbar = eos->GetZbar(index, rho, Te);
   double nu_ei = ClassicalMeanStoppingPower::Eval(T, ip, rho);
   // In classical approach nu_ei = Zbar*nu_ee is assumed.
   double nu_ee = nu_ei / Zbar;
   // AWBS correction based on comparison of diffusive asymptotic of AWBS 
   // and SH, and a resulting dependence on Zbar.
   //corrAWBS = (688.9*Zbar + 114.4) / (Zbar*Zbar + 1038.0*Zbar + 474.1);
   //corrAWBS = 1.0;

   return GetCorrAWBS(Zbar) * nu_ee;
}

void LorentzEfield::Eval(Vector &V, ElementTransformation &T,
                           const IntegrationPoint &ip)
{
   double rho = rho_gf.GetValue(T.ElementNo, ip);
   double Te = Te_gf.GetValue(T.ElementNo, ip);
   double index = material_pcf->Eval(T, ip);
   double Zbar = eos->GetZbar(index, rho, Te);
   //double Te = max(1e-10, Te_gf.GetValue(T.ElementNo, ip));
   // Set the scaled velocity to correspond to the local thermal velocity.
   double vTe = eos->GetvTe(Te);
   // Compute the temperature and density length scales.
   Vector grad_Te;
   Te_gf.GetGradient(T, grad_Te);
   Vector grad_rho;
   rho_gf.GetGradient(T, grad_rho);
   // Return the Lorentz quasi-neutral (zero current) Efield.
   V = grad_Te;
   // Efield is scaled as Escaled = q/me*E, which is also used in Lorentz force.
   // Original Lorentz zero-current condition.
   double xi = 2.5;
   // SH E field correction.
   //double xi = 1.0 + 1.5 * (Zbar + 0.477) / (Zbar + 2.15); 
   V *= S0 * pow(vTe, 2.0) * xi / Te;
/*
   V = grad_Te;
   V *= 2.5 / Te;
   V *= rho;
   V += grad_rho;
   V *= 1.0 / rho;
   // S0 acts as electric field scaling, if necessary.
   V *= S0 * pow(vTe, 2.0);
*/
}

double LorentzEfield::Eval(ElementTransformation &T,
                           const IntegrationPoint &ip)
{
   // Return the Lorentz quasi-neutral (zero current) Efield.
   Vector V;
   Eval(V, T, ip); 

   return V.Norml2();
}

double P1a0KineticCoefficient::Eval(ElementTransformation &T,
                                    const IntegrationPoint &ip)
{
   // "velocity_real" must be set previously.
   mspei_pcf->SetVelocityReal(velocity_real);
   mspee_pcf->SetVelocityReal(velocity_real);
   double nu_ei =  mspei_pcf->Eval(T, ip);
   double nu_ee =  mspee_pcf->Eval(T, ip);

   double dFMdv = dFM->GetValue(T.ElementNo, ip);
   double dF0dv = dF0->GetValue(T.ElementNo, ip); 

   double nuEe_scale, nuEt_scale, Efield_scale;
   Efield_pcf->GetEscales(T, ip, velocity_real, nu_ee,
                          nuEe_scale, nuEt_scale, Efield_scale);

   double nu_t = nu_ei + (1.0 + nuEt_scale) * nu_ee;
   //double nu_t = nu_ei + nu_ee + nu_E;

   // Scattering on ions and electrons.
   //return velocity_real / 3.0 / nu_t * dF0dv * velocity_real * velocity_real;
   return velocity_real / 3.0 / nu_t * (dFMdv + Efield_scale * dF0dv) 
          * velocity_real * velocity_real;
}

void P1b0KineticCoefficient::Eval(Vector &V, ElementTransformation &T,
                                  const IntegrationPoint &ip)
{
   // "velocity_real" must be set previously.
   mspei_pcf->SetVelocityReal(velocity_real);
   mspee_pcf->SetVelocityReal(velocity_real);
   double nu_ei =  mspei_pcf->Eval(T, ip);
   double nu_ee =  mspee_pcf->Eval(T, ip);

   double nu_E, nuEe_scale, nuEt_scale, Efield_scale;
   Efield_pcf->GetEscales(T, ip, velocity_real, nu_ee,
                          nuEe_scale, nuEt_scale, Efield_scale);
   // Provide appropriate isotropic Efield effect.
   nu_E = nu_ee * nuEe_scale;

   double nu_t = nu_ei + (1.0 + nuEt_scale) * nu_ee;
   //double nu_t = nu_ei + nu_ee + nu_E;

   T.SetIntPoint(&ip);
   F0->GetGradient(T, V); 
   // Scattering on ions and electrons.
   V *= 1.0 / 3.0 / nu_t;
   
   for (int d = 0; d < vdim; d++)
   { 
      // Scattering on ions and electrons.
      V(d) -= (nu_ee + nu_E) / nu_t * dF1->GetValue(T.ElementNo, ip, d);
      //V(d) -= nu_ee / nu_t * dF1->GetValue(T.ElementNo, ip, d);
   }
   
   V *= pow(velocity_real, 4.0);
}

void P1b1KineticCoefficient::Eval(Vector &V, ElementTransformation &T,
                                  const IntegrationPoint &ip)
{
   // "velocity_real" must be set previously.
   mspei_pcf->SetVelocityReal(velocity_real);
   mspee_pcf->SetVelocityReal(velocity_real);
   double nu_ei =  mspei_pcf->Eval(T, ip);
   double nu_ee =  mspee_pcf->Eval(T, ip);

   double nu_t = nu_ei + nu_ee;

   for (int d = 0; d < vdim; d++)
   { 
      V(d) = F1->GetValue(T.ElementNo, ip, d);
   }
   
   // Scattering on ions and electrons.
   V *= pow(velocity_real, 3.0) / nu_t;
}

void OhmCurrentCoefficient::Eval(Vector &V, ElementTransformation &T,
                                 const IntegrationPoint &ip)
{
   Vector Ex(vdim), Bx(vdim);
   Efield_pcf->Eval(Ex, T, ip);
   Bfield_pcf->Eval(Bx, T, ip);

   double a0x = a0_pgf->GetValue(T.ElementNo, ip);
   for (int d = 0; d < vdim; d++)
   { 
	  double b0x_d = b0_pgf->GetValue(T.ElementNo, ip, d);
      double b1x_d = b1_pgf->GetValue(T.ElementNo, ip, d);
	  V(d) = - b0x_d - a0x * Ex(d);
   }
}

void OhmEfieldCoefficient::Eval(Vector &V, ElementTransformation &T,
                                const IntegrationPoint &ip)
{
   Vector Bx(vdim);
   Bfield_pcf->Eval(Bx, T, ip);

   double a0x = a0_pgf->GetValue(T.ElementNo, ip);
   for (int d = 0; d < vdim; d++)
   { 
	  double jCx_d = jC_pgf->GetValue(T.ElementNo, ip, d);
	  double b0x_d = b0_pgf->GetValue(T.ElementNo, ip, d);
      double b1x_d = b1_pgf->GetValue(T.ElementNo, ip, d);
	  V(d) = - b0x_d / a0x - alpha * jCx_d / a0x;
   }
}

double MaxwellBoltzmannCoefficient::Eval(ElementTransformation &T,
                                         const IntegrationPoint &ip, double rho)
{
   // "velocity_real" must be set previously.
   double pi = 3.14159265359;
   double Te = max(1e-10, Te_gf.GetValue(T.ElementNo, ip));
   double vTe = eos->GetvTe(Te);
   double index = material_pcf->Eval(T, ip);
   double ne = eos->GetElectronDensity(index, rho, Te);

   // Maxwell-Boltzmann distribution fM = ne*vT^3*(2/pi)^1.5*exp(-v^2/2/vT^2)
   double fM = 4.0 * pi *
               ne / pow(vTe, 3.0) / pow(2.0 * pi, 1.5) *
               exp(- pow(velocity_real, 2.0) / 2.0 / pow(vTe, 2.0));

   return S0 * fM;
}

double MaxwellBoltzmannCoefficient::Eval(ElementTransformation &T, 
                                         const IntegrationPoint &ip)
{
   double rho = rho_gf.GetValue(T.ElementNo, ip);

   return Eval(T, ip, rho);
}

double AWBSdfMdv::Eval(ElementTransformation &T,
                          const IntegrationPoint &ip, double rho)
{
   // "velocity_real" must be set previously.
   double pi = 3.14159265359;
   double Te = max(1e-10, Te_gf.GetValue(T.ElementNo, ip));
   double vTe = eos->GetvTe(Te);
   //double index = material_pcf->Eval(T, ip);
   //double ne = eos->GetElectronDensity(index, rho, Te);

   // Maxwell-Boltzmann distribution fM = ne*vT^3*(2/pi)^1.5*exp(-v^2/2/vT^2)
   //double fM = 4.0 * pi *
   //            ne / pow(vTe, 3.0) / pow(2.0 * pi, 1.5) *
   //            exp(- pow(velocity_real, 2.0) / 2.0 / pow(vTe, 2.0));
   MaxBo_cf.SetVelocityReal(velocity_real); 
   double fM = MaxBo_cf.Eval(T, ip, rho);
   double dfMdv = - velocity_real / pow(vTe, 2.0) * fM;

   // S0 acts as electron source scaling (scaling electron density), 
   // if necessary.
   return S0 * dfMdv;
}

double AWBSdfMdv::Eval(ElementTransformation &T, const IntegrationPoint &ip)
{
   double rho = rho_gf.GetValue(T.ElementNo, ip);

   return Eval(T, ip, rho);
}

// Diagnostics.
double CorrEfieldCoefficient::Eval(ElementTransformation &T,
                                   const IntegrationPoint &ip)
{
   // Multiple of the local thermal velocity, 
   // i.e. the maximum velocity considered is N_x * vTe.
   int N_x = 99;
   int Nd = 1000;
   double vTe = mspei_pcf->GetvTe(T, ip);
   double dv = N_x * vTe / Nd;

   Vector Efield(vdim);
   Efield = 1.0;
   Efield_pcf->Eval(Efield, T, ip);
   double Efield_norm = Efield.Norml2();

   // Start with a low velocity compared to vTe.
   double velocity_limit = Nd * dv;
   for (int i = Nd; i > 1; i--)
   {
      double velocity_real = i*dv;
      mspee_pcf->SetVelocityReal(velocity_real);
      double nu_ee =  mspee_pcf->Eval(T, ip);

      double Escale = 1.0;
      double corrE = Escale * nu_ee * velocity_real / Efield_norm;
      double Efield_portion = 0.5;
	  if (corrE < Efield_portion)
      //if (scale * nu_ee * velocity_real < Efield_norm) 
      {
         velocity_limit = velocity_real;
         //break;
      }
   }

   //cout << "velocity_limit/vTe: " << velocity_limit / vTe << endl;

   // Return velocity limit scaled to local thermal velocity.
   return velocity_limit / vTe;
}


} // namespace nth

} // namespace mfem

#endif // MFEM_USE_MPI

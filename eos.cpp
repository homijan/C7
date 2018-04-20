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

#include "eos.hpp"

#ifdef MFEM_USE_MPI

using namespace std;

namespace mfem
{

namespace nth
{

// Objects for input profile data reading.
InputProfile *inElectronTemp = NULL, *inElectronDens = NULL, *inZbar = NULL;

// Constant profile contructor.
InputProfile::InputProfile(double _const_data)
{
   // Make sure that the vector data are empty.
   x_data.resize(0);
   data.resize(0);
   // Store the constant value of the profile.
   const_data = _const_data;
}

InputProfile::InputProfile(std::string filename)
{
   x_data.resize(0);
   data.resize(0);
   double x, value;
   //std::ifstream ifs("VFPdata/temperature.dat", std::ifstream::in);
   std::ifstream ifs(filename.c_str(), std::ifstream::in);
   while(!ifs.eof())
   {
      ifs >> x;
      ifs >> value;
      // Store x coordinate in [cm].
      x_data.push_back(x);
      //x_data.push_back(x * 1e-4);
      // Store temperature in [eV].
      data.push_back(value); 
      //data.push_back(value * 1e3);
   }
}

double InputProfile::GetValue(double x)
{
   // Constant profile.
   if (data.size() == 0) { return const_data; }
   // Data input profile.
   int i0, i1, ic;
   // Half interval method.
   i0 = 0;
   i1 = data.size() - 1;
   // Keep x in the range.
   if (x <= x_data[i0]) { return data[i0]; }
   if (x >= x_data[i1]) { return data[i1]; }
   // Find the right interval.
   while (i1 - i0 > 1)
   { 
      ic = (i1 + i0) / 2;
	  if ((x - x_data[ic])*(x - x_data[i0]) <= 0) { i1 = ic; }
      else { i0 = ic; }
	  //cout << "i0, i1: " << i0 << ", " << i1 << endl;
   }  
   return data[i0] + (data[i1] - data[i0]) / (x_data[i1] - x_data[i0])
          * (x - x_data[i0]);
   
   /*
   int index = 0;
   if (x <= x_data[0]) { return data[0]; }
   if (x >= x_data[x_data.size()-1]) { return data[data.size()-1]; }
   for (int i = 0; i < x_data.size() - 1; i++)
   {
      if (x > x_data[i] && x <= x_data[i+1])
      {
         return data[i] + (data[i+1] - data[i]) / (x_data[i+1] - x_data[i])
                 * (x - x_data[i]);
      }
   }
   */
}

double ExtDataEOS::GetX(double rho, double Te)
{
   double xTe = 0.0;
   vector<double> &x = inElectronTemp->x_data;
   vector<double> &Td = inElectronTemp->data;
   int s = inElectronTemp->data.size();
   int i0, i1, ic;
   double Te0, Te1, x0, x1;
   // Half interval method.
   i0 = 0;
   i1 = s-1;
   // Keep x in the range.
   if (Td[i0] > Td[i1] && Te > Td[i0]) { return x[i0]; }
   if (Td[i0] > Td[i1] && Te < Td[i1]) { return x[i1]; }
   if (Td[i0] < Td[i1] && Te < Td[i0]) { return x[i0]; }
   if (Td[i0] < Td[i1] && Te > Td[i1]) { return x[i1]; }
   // Find the right interval.
   while (i1 - i0 > 1)
   { 
      ic = (i1 + i0) / 2;
	  if ((Te - Td[ic])*(Te - Td[i0]) <= 0) { i1 = ic; }
      else { i0 = ic; }
	  //cout << "i0, i1: " << i0 << ", " << i1 << endl;
   }
   // Linear data approximation.
   Te0 = Td[i0];
   Te1 = Td[i1];
   x0 = x[i0];
   x1 = x[i1];
   xTe = (Te - Te0) / (Te1 - Te0) * (x1 - x0) + x0;
   //cout << "Te, Te0, Te1, xTe: " << Te << ", " << Te0 << ", " << Te1 << ", " 
   //     << xTe << endl;
   return xTe;
   
   /*
   for (int i = 0; i < s - 1; i++)
   {
      Te0 = Td[i];
      Te1 = Td[i+1];
      x0 = x[i];
      x1 = x[i+1];
      if ((Te >= Te0 && Te <= Te1) || (Te <= Te0 && Te >= Te1))
      {
         xTe = (Te - Te0) / (Te1 - Te0) * (x1 - x0) + x0;
      }
   }

   // TMP approach considering monotone Te data.
   return xTe;
   */
}

double ExtDataEOS::GetElectronDensity(double index, double rho, double Te)
{
   // Use constant value.
   if (inElectronDens == NULL) { return const_ne; }
   // Use data value based on x coordinate.
   double x = GetX(rho, Te);
   return inElectronDens->GetValue(x);
}

double ExtDataEOS::GetZbar(double index, double rho, double Te)
{
   // Use constant value.
   if (inZbar == NULL) { return const_Zbar; } 
   // Use data value based on x coordinate.
   double x = GetX(rho, Te);
   return inZbar->GetValue(x);
}

double ExtDataEOS::GetIonDensity(double index, double rho, double Te)
{
   return GetElectronDensity(index, rho, Te) / GetZbar(index, rho, Te);
}

} // namespace nth

} // namespace mfem

#endif // MFEM_USE_MPI

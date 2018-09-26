import numpy as np
## Math formulas.
from math import pi
from math import exp
## Interpolation.
from scipy.interpolate import splev, splrep
## Brent solver for corr factors.
from scipy import optimize
## delta corr fit.
from scipy.optimize import curve_fit
## Graphics.
import matplotlib.pyplot as plt
import matplotlib

## Global setting of plotting.
font = {'family' : 'Sans',
        #'weight' : 'bold',
        'size'   : 19}
figure = {'figsize' : '10.5, 6.5'} # a tuple in inches
matplotlib.rc('font', **font)
matplotlib.rc('figure', **figure)

## Fundamental physical constants in cgs. 
kB = 1.6022e-12
me = 9.1094e-28
qe = 4.8032e-10

def vTh(T): 
    return (kB*T/me)**0.5
def fM(n, T, v):
    return n/(vTh(T)**3.0*(2.0*pi)**1.5)*exp(-v**2.0/2.0/vTh(T)**2.0)
def dfMdz(n, T, dndz, dTdz, v):
    dfMdn = fM(n, T, v) / n
    dfMdT = (0.5*v**2.0 / vTh(T)**2.0 - 3.0 / 2.0) / T * fM(n, T, v)
    return dfMdn * dndz + dfMdT * dTdz

## Cross section in cgs.
sigma = 8.1027575e17 ## Matching the SH diffusive flux.
## Coulomb logarithm.
coulLog = 7.09
coulLogLv = 7.45 


## Some default values.
ne = 1.0e23
Te = 100.0
dTedz = -1e0
dnedz = 0.0
Zbar = 1.0
ZbarLv = 116.0

## Full electron-electron cross section.
Gamma_ee = sigma * coulLog
GammaLv_ee = sigma * coulLogLv
#Gamma_ee = 1.0

lambda_ei = vTh(Te)**4.0 / Gamma_ee / ne / (Zbar + 1.0)
Kn_ei = lambda_ei * dTedz / Te
print "Kn_ei: ", Kn_ei 

import argparse
## Create parser object.
parser = argparse.ArgumentParser(description='Compare the diffusive asymptotic to C7 computation.')
## Define input arguments.
parser.add_argument("-s", "--sigma", help="Sigma for electro-ion cross-section.", type=float)
parser.add_argument("-Z", "--Zbar", help="Zbar for the ion charge.", type=float)
parser.add_argument("-Np", "--Nproc", help="Number of processors used to compute the data.", type=int)
## A no value argument solution.
parser.add_argument("-ps", "--pltshow", action='store_true', help="Plot show() by adding -ps/--pltshow argument.")
parser.add_argument("-lF1", "--labelFluxExt1", help="Use -lF1/--labelFluxExt1 to use and label VFPdata/flux1.dat.")


## Parse arguments.
args = parser.parse_args()
if (args.sigma):
   sigma = args.sigma
if (args.Zbar):
   Zbar = args.Zbar
###############################################################################
########### Analysis of diffusive asymptotic of AWBS model #################### 
###############################################################################
############################################################################### 

###############################################################################
#### FP equation diffusive regime #############################################
def RosenbluthPotentialsF0(vs, f0s):
    ## Find number of cells and use equidistant velocity step.
    N = len(vs)
    dv = (max(vs) - min(vs)) / (N - 1)
    ## Define Rosenbluth potential arrays.
    I0 = np.zeros(N)
    I2 = np.zeros(N)
    Jm1 = np.zeros(N) 
    ## Integrate Rosenbluth potentials with zero starting values.
    ## Ascending potentials.
    I0[0] = 0.0
    I2[0] = 0.0 
    for i in range(1, N):
        #print i
        I0[i] = I0[i-1] + (f0s[i-1]*vs[i-1]**2.0 + f0s[i]*vs[i]**2.0) * dv / 2.0
        I2[i] = I2[i-1] + (f0s[i-1]*vs[i-1]**4.0 + f0s[i]*vs[i]**4.0) * dv / 2.0
    for i in range(0, N):
        I2[i] = I2[i] / vs[i]**2.0
    ## Descending potentials.
    Jm1[N-1] = 0.0
    for i in range(N-2, -1, -1):
        #print i
        Jm1[i] = Jm1[i+1] + (f0s[i+1]*vs[i+1] + f0s[i]*vs[i]) * dv / 2.0
    for i in range(N-2, -1, -1):
        Jm1[i] = Jm1[i] * vs[i]
    ## The raight scale.
    I0 = 4.0 * pi * I0
    I2 = 4.0 * pi * I2
    Jm1 = 4.0 * pi * Jm1        
    return I0, I2, Jm1 

def RosenbluthPotentialsF1(vs, f1s):
    ## Find number of cells and use equidistant velocity step.
    N = len(vs)
    dv = (max(vs) - min(vs)) / (N - 1)
    ## Define Rosenbluth potential arrays.
    I1 = np.zeros(N)
    I3 = np.zeros(N)
    Jm2 = np.zeros(N) 
    ## Integrate Rosenbluth potentials with zero starting values.
    ## Ascending potentials.
    I1[0] = 0.0
    I3[0] = 0.0 
    for i in range(1, N):
        #print i
        I1[i] = I1[i-1] + (f1s[i-1]*vs[i-1]**3.0 + f1s[i]*vs[i]**3.0) * dv / 2.0
        I3[i] = I3[i-1] + (f1s[i-1]*vs[i-1]**5.0 + f1s[i]*vs[i]**5.0) * dv / 2.0
    for i in range(0, N):
        I1[i] = I1[i] / vs[i]
        I3[i] = I3[i] / vs[i]**3.0
    ## Descending potentials.
    Jm2[N-1] = 0.0
    for i in range(N-2, -1, -1):
        #print i
        Jm2[i] = Jm2[i+1] + (f1s[i+1] + f1s[i]) * dv / 2.0
    for i in range(N-2, -1, -1):
        Jm2[i] = Jm2[i] * vs[i]**2.0
    ## The raight scale.
    I1 = 4.0 * pi * I1
    I3 = 4.0 * pi * I3
    Jm2 = 4.0 * pi * Jm2        
    return I1, I3, Jm2

def FPcoefficientsF0(vs, f0s):
    ## Compute Rosenbluth potentials using f0.
    I0, I2, Jm1 = RosenbluthPotentialsF0(vs, f0s)
    ## Fill the collision operator C_ee coefficients using f0.
    C_d2fdv2 = (I2 + Jm1) / 3.0 / vs
    C_dfdv = (3.0 * I0 - I2 + 2.0 * Jm1) / 3.0 / vs**2.0
    Ce_f = 8.0*pi*fMs - (3.0 * I0 - I2 + 2.0 * Jm1) / 3.0 / vs**3.0
    ## Fill the collision operator C_ei coefficient using delta_i.
    Ci_f = - Zbar * ne / vs**3.0
    return C_d2fdv2, C_dfdv, Ce_f, Ci_f

def d2fdv2_dfdv(vs, fs):
    ## Compute the first and second derivatives from data interpolation.
    smooth = 0 # lower less smoothing
    ## Find a spline for the f0 data.
    ftck = splrep(vs, fs, s=smooth)
    dfdv = splev(vs, ftck, der=1)
    d2fdv2 = splev(vs, ftck, der=2)
    return d2fdv2, dfdv

def FPtermsF0F1(vs, f0s, f1s, df0dzs, Ez, Gamma):
    ## Compute Rosenbluth potentials using f1.
    I1, I3, Jm2 = RosenbluthPotentialsF1(vs, f1s)
    ## Compute the first and second derivative of f0.
    d2f0dv2, df0dv = d2fdv2_dfdv(vs, f0s)
    ## Compute the explicit term of FP using f0 and f1.
    d2f0dv2If1_df0dvIf1 = - (1.0 / 5.0 / vs * (I3 + Jm2) * d2f0dv2 + 1.0 / 15.0 / vs**2.0 * (5.0 * I1 - 3.0 * I3 + 2.0 * Jm2) * df0dv)
    # Advection part of FP using f0.
    df0dz = 1.0 / Gamma * vs * df0dzs
    Edf0dv = 1.0 / Gamma * Ez * df0dv
    return d2f0dv2If1_df0dvIf1, df0dz, Edf0dv

def FPtermsFMF1(vs, fMs, f1s, df0dzs, Ez, T, Gamma):
    ## Compute Rosenbluth potentials using f1.
    I1, I3, Jm2 = RosenbluthPotentialsF1(vs, f1s) 
    ## Compute the explicit term of FP using fM and f1.
    d2fMdv2If1_dfMdvIf1 = - fMs * 1.0 / 15.0 / vs / vTh(T)**2.0 * (3.0 * vs**2.0 / vTh(T)**2.0 * (I3 + Jm2) - 5.0 * (I1 + Jm2))
    # Advection part of FP using fM.
    dfMdz = 1.0 / Gamma * vs * dfMdzs
    EdfMdv = - 1.0 / Gamma * vs * Ez / vTh(T)**2.0 * fMs
    return d2fMdv2If1_dfMdvIf1, dfMdz, EdfMdv

def AllFokkerPlanckEquationTerms(vs, f0s, f1s, df0dzs, Ez, Gamma):
    ## The Fokker-Planck equation for electrons in diffusive regime
    ## can be written as
    ##
    ## C2(f0) d2f1dv2 + C1(f0) * df1dv + (C0e(f0) + Ci) * f1 = 
    ## B2(f1) * d2f0dv2 + B1(f1) * df0dv + v * df0dz + Ez * df0dv
    ##
    ## Compute the coefficients of f1 terms.
    C2, C1, C0e, C0i = FPcoefficientsF0(vs, f0s)
    ## Compute the first and second derivative of f1.
    d2f1dv2s, df1dvs = d2fdv2_dfdv(vs, f1s)
    ## Assemble the f1 terms.
    C2d2f1dv2 = C2 * d2f1dv2s
    C1df1dv = C1 * df1dvs
    C0f1 = (C0e + C0i) * f1s
    ## Compute the f0 terms.
    B2d2f0dv2If1_B1df0dvIf1, df0dz, Edf0dv = FPtermsF0F1(vs, f0s, f1s, df0dzs, Ez, Gamma)
    ## Assemble the f0 terms.
    B2d2f0dv2_B1df0dv_df0dz_Edf0dv = B2d2f0dv2If1_B1df0dvIf1 + df0dz + Edf0dv
    return C2d2f1dv2, C1df1dv, C0f1, B2d2f0dv2_B1df0dv_df0dz_Edf0dv

def bdf_cf_equals_d(vs, b, c, d, f0):
    ## Solving equation b*dfdv + c*f = d
    ## Create discretization of f.
    ## Find number of cells and use equidistant velocity step.
    N = len(vs)
    dv =  (max(vs) - min(vs)) / (N - 1)
    ## Define the unknowns.
    f = np.zeros(N)
    ## Set initial conditions for infinity->0 integration.
    f[N-1] = f0
    # Perform backward Euler integration.   
    for i in range(N-2, -1, -1):
        bb = b[i+1]
        cc = c[i+1]
        dd = d[i+1]
        #bb = (b[i] + b[i+1]) / 2.0
        #cc = (c[i] + c[i+1]) / 2.0
        #dd = (d[i] + d[i+1]) / 2.0
        f[i] = (dd + bb / dv * f[i+1]) / (cc + bb / dv) 
    return f

def ad2f_bdf_cf_equals_d(vs, a, b, c, d, f0, df0):
    ## Solving equation a*dfdv2 + b*dfdv + c*f = d
    ## omega = sqrt(a / c), theta = b / 2 / sqrt(a * c)
    ## Create discretization of f.
    ## Find number of cells and use equidistant velocity step.
    N = len(vs)
    dv =  (max(vs) - min(vs)) / (N - 1)
    #omega = (a / c)**0.5
    #theta = b / 2.0 / (a * c)**0.5
    #print "omega, theta, omega/dv: ", omega[0], theta[0], omega[0] / dv
    ## Define the unknowns.
    f = np.zeros(N)
    df = np.zeros(N)
    ## Set initial conditions for infinity->0 integration.
    f[N-1] = f0
    df[N-1] = df0
    # Perform backward Euler integration.   
    for i in range(N-2, -1, -1):
        aa = (a[i] + a[i+1]) / 2.0
        bb = (b[i] + b[i+1]) / 2.0
        cc = (c[i] + c[i+1]) / 2.0
        dd = (d[i] + d[i+1]) / 2.0
        df[i] = (dv * dd + aa * df[i+1] - dv * cc * f[i+1]) / (dv**2.0 * cc + aa + dv * bb)
        f[i] = (dv * dd + aa * df[i+1] - (aa + dv * bb) * df[i]) / dv / cc 
    return f, df
#### FP equation diffusive regime #############################################
###############################################################################

###############################################################################
########### AWBS diffusive asymptotic ######################################### 
#def AWBS_distribution(v, f0, ne, Te, gradTe, Z, E, Gamma_ee, corr):
#    N = len(v)
#    f1 = np.zeros(N) 
#    j1 = np.zeros(N)
#    q1 = np.zeros(N)
#    f1[N-1] = f0
#    for i in range(N-1, 0, -1):
#        dv = v[i] - v[i-1]
#        vp = v[i-1]
#        mfpei = vp**4.0 / Gamma_ee / ne / Z
#        rhs = corr * Z * mfpei / (corr * Z + 1.0) * ((vp**2.0 / 2.0 / vTh(Te)**2.0 - 1.5) * gradTe / Te - E / vTh(Te)**2.0)
#        rhs = rhs * fM(ne, Te, vp)
#        f1[i-1] = (f1[i] * vp / (corr * Z + 1.0) / dv - rhs) /  (1.0 + vp / (corr * Z + 1.0) / dv) # ee iso 
#        j1[i-1] = f1[i-1] * vp**3.0
#        q1[i-1] = f1[i-1] * vp**5.0
#    return f1, j1, q1
def _AWBS_distribution(v, f0, ne, Te, gradTe, Z, E, Gamma_ee, corr):
    N = len(v)
    f1 = np.zeros(N) 
    j1 = np.zeros(N)
    q1 = np.zeros(N)
    f1[N-1] = f0
    for i in range(N-1, 0, -1):
        dv = v[i] - v[i-1]
        vp = v[i-1]
        mfpe = vp**4.0 / Gamma_ee / ne
        rhs = corr * mfpe / vp * ((vp**2.0 / 2.0 / vTh(Te)**2.0 - 1.5) * gradTe / Te - E / vTh(Te)**2.0)
        rhs = rhs * fM(ne, Te, vp)
        f1[i-1] = (- f1[i] / dv + rhs) /  (- (corr * Z + 1.) / vp - 1. / dv) 
        j1[i-1] = f1[i-1] * vp**3.0
        q1[i-1] = f1[i-1] * vp**5.0
    return f1, j1, q1
def HighVelocity_distribution(v, f0, ne, Te, gradTe, Z, E, Gamma_ee, corr):
    N = len(v)
    f1 = np.zeros(N) 
    j1 = np.zeros(N)
    q1 = np.zeros(N)
    f1[N-1] = f0
    for i in range(N-1, 0, -1):
        dv = v[i] - v[i-1]
        vp = v[i-1]
        mfpe = vp**4.0 / Gamma_ee / ne
        rhs = corr * mfpe / vp * ((vp**2.0 / 2.0 / vTh(Te)**2.0 - 1.5) * gradTe / Te - E / vTh(Te)**2.0)
        rhs = rhs * fM(ne, Te, vp)
        f1[i-1] = (- f1[i] / dv + rhs) /  (- (Z + 1. / 2.) / vp - 1. / dv) 
        j1[i-1] = f1[i-1] * vp**3.0
        q1[i-1] = f1[i-1] * vp**5.0
    return f1, j1, q1

## Correct AWBS diffusive limit distribution calculation! 
## Using appropriate explicit values of corr_nu and corr_E.
def AWBS_distribution(v, ne, Te, gradTe, Z, E, Gamma_ee, corr_nu=0.5, corr_E=1.0):
    # corr stands for nue^* = corr * nue, i.e. mfpe^* = mfpe / corr. 
    # corr * v / mfpe * df1dv - (Z + corr) / mfpe * f1 = 
    # dfMdz + qe * E / me / v * dfMdv
    N = len(v)
    f1 = np.zeros(N) 
    j1 = np.zeros(N)
    q1 = np.zeros(N)
    # Integration starts from zero particle density for velocity -> infinity.
    f0 = 0.0 
    f1[N-1] = f0
    for i in range(N-1, 0, -1):
        dv = v[i] - v[i-1]
        vp = v[i-1]  
        mfpe = vp**4.0 / Gamma_ee / ne
        rhs = mfpe / corr_nu / vp * ((vp**2.0 / 2.0 / vTh(Te)**2.0 - 1.5) * gradTe / Te - corr_E * E / vTh(Te)**2.0)
        rhs = rhs * fM(ne, Te, vp)
        f1[i-1] = (- f1[i] / dv + rhs) /  (- (Z + corr_nu) / corr_nu / vp - 1. / dv) 
        j1[i-1] = f1[i-1] * vp**3.0
        q1[i-1] = f1[i-1] * vp**5.0 
    return f1, j1, q1

## AWBS distribution calculation in order to seek for j=0 based on corr_E.
def corr_E_AWBS_distribution(corr_E, v, ne, Te, gradTe, Z, E, Gamma_ee, corr_nu):
    ## Compute the q1 AWBS distribution on velocities v. 
    f1, j1, q1 = AWBS_distribution(v, ne, Te, gradTe, Z, E, Gamma_ee, corr_nu, corr_E)
    ## Integrate heat flux along equidistant velocity.
    dv = v[1] - v[0]
    j = qe * 4.0 / 3.0 * np.pi * sum(j1) * dv    
    # Zero current condition.
    return j

## AWBS distribution calculation in order to seek for q_AWBS = q_SH 
## based on corr_nu. The zero current condition j=0 holds.
def corr_nu_AWBS_distribution(corr_nu, v, ne, Te, gradTe, Z, E, Gamma_ee, q_SH):
    #corr_E = 1.0
    corr_E =  optimize.brentq(corr_E_AWBS_distribution, 0.1, 10.0, args=(v, ne, Te, gradTe, Z, E, Gamma_ee, corr_nu))
    #print "corr_E:", corr_E
    ## Compute the q1 AWBS distribution on velocities v. 
    f1, j1, q1 = AWBS_distribution(v, ne, Te, gradTe, Z, E, Gamma_ee, corr_nu, corr_E)
    ## Integrate heat flux along equidistant velocity.
    dv = v[1] - v[0]
    q_AWBS = me / 2.0 * 4.0 / 3.0 * np.pi * sum(q1) * dv    
    return abs(q_AWBS) - abs(q_SH)

## Usefull overall output from AWBS distribution calculation 
## with respect to original SH1953 calculations.
def FindAWBSdependenceOnZ(N, ne, Te, gradTe, Gamma_ee): 
    ## Number of cells in velocity magnitude.
    #N = 10000
    ## Multiple of thermal velocity for min(v).
    min_x_vTh = 0.001
    ## Multiple of thermal velocity for max(v).
    max_x_vTh = 7.0
    ## Define velocity magnitude discretization.
    v = np.linspace(min_x_vTh * vTh(Te), max_x_vTh * vTh(Te), N)
    dv = v[1] - v[0]
    ## Lorentz gas E field valid for Z>>1.
    E_Lorentz = vTh(Te)**2.0 * (dnedz / ne +  2.5 * dTedz / Te)

    ## Find AWBS correction to nu and E for the following Z.
    #Zs = [1.0, 2.0, 4.0, 16.0, 20.0, 50.0, 116.0, 1000.0]
    Zs = [1.0, 2.0, 4.0, 10.0, 16.0]
    corrs = []
    for Z in Zs:
        print "Z:", Z
        ## The case Z = 1 provides good distribution function in SH paper.
        ## In other cases the distribution function is not so convincing.
        # Get the original SH heat flux. 
        #f1SH1953, j1SH1953, q1SH1953 = SH_distribution(v, ne, Te, gradTe, Z, Gamma_ee)
        #q_SH1953 = me / 2.0 * 4.0 / 3.0 * np.pi * sum(q1SH1953) * dv
        
        #Gamma_ee = sigma * coulLog
        coulLog = Gamma_ee / sigma  
        SHcorr = SH_corr(Z)
        ## SH computation was FP based up to Z = 4.
        if (Z > 4):
            SHcorr = (Z + 0.24)/(Z + 4.2)  
        q_SH = - SHcorr * 1.31e10 / coulLog / Z * Te**2.5 * gradTe
        #print "q_SH1953:", q_SH1953
        #print "q_SH:", q_SH 
        ## Find appropriate AWBS corr_nu matching the SH heat flux along j=0.
        corr_nu =  optimize.brentq(corr_nu_AWBS_distribution, 0.01, 1000.0, args=(v, ne, Te, gradTe, Z, E_Lorentz, Gamma_ee, q_SH))
        corrs.append(corr_nu - 0.5)
        print "corr_nu:", corr_nu 
       
        ## Obtain the appropriate corr_E.
        corr_E =  optimize.brentq(corr_E_AWBS_distribution, 0.1, 10.0, args=(v, ne, Te, gradTe, Z, E_Lorentz, Gamma_ee, corr_nu))
        print "corr_E:", corr_E

        ## Double-check the corr_nu and corr_E results.
        f1, j1, q1 = AWBS_distribution(v, ne, Te, gradTe, Z, E_Lorentz, Gamma_ee, corr_nu, corr_E)
        j_AWBS = qe * 4.0 / 3.0 * np.pi * sum(j1) * dv
        q_AWBS = me / 2.0 * 4.0 / 3.0 * np.pi * sum(q1) * dv
        print "j_AWBS:", j_AWBS
        print "(q_AWBS - q_SH) / q_SH:", abs((q_AWBS - q_SH) / q_SH)

        ## General case of using corr_nu = 0.5 always.
        dcorr = corr_nu - 0.5
        corr_nu = 0.5
        print "The case of always using corr_nu = 0.5!"
        print "corr_nu deviation from 0.5 =", dcorr
        ## Obtain the appropriate corr_E.
        corr_E =  optimize.brentq(corr_E_AWBS_distribution, 0.1, 10.0, args=(v, ne, Te, gradTe, Z, E_Lorentz, Gamma_ee, corr_nu))
        print "corr_E:", corr_E

        ## Double-check the corr_nu and corr_E results.
        f1, j1, q1 = AWBS_distribution(v, ne, Te, gradTe, Z, E_Lorentz, Gamma_ee, corr_nu, corr_E)
        j_AWBS = qe * 4.0 / 3.0 * np.pi * sum(j1) * dv
        q_AWBS = me / 2.0 * 4.0 / 3.0 * np.pi * sum(q1) * dv
        print "j_AWBS:", j_AWBS
        print "(q_AWBS - q_SH) / q_SH:", abs((q_AWBS - q_SH) / q_SH) 

    ## Find an analytic fit to corr_nu.
    Z = np.linspace(1, 100, 1000)  
    popt, pcov = curve_fit(rational2_2, Zs, corrs)
    u0 = popt[1]
    u1 = popt[0]
    l0 = popt[3]
    l1 = popt[2]   
    print "fit = (u0 + u1 * Z) / (l0 + l1 * Z), u0 ,u1, l0, l1:", u0, u1, l0, l1

    plt.plot(Zs, corrs, label='corr - 0.5')
    #plt.plot(Z, (u0 + u1 * Z) / (l0 + l1 * Z), label='delta fit')
    plt.plot(Z, (-1.11 + 0.59 * Z) / (5.15 + 8.37 * Z), label='fit') 

    #plt.plot(Z, rational2_2(Z, *popt), label='fit')

    plt.legend()
    plt.show()
    return corr_nu

def rational(x, p, q):
    """
    The general rational function description.
    p is a list with the polynomial coefficients in the numerator
    q is a list with the polynomial coefficients (except the first one)
    in the denominator
    The zeroth order coefficient of the denominator polynomial is fixed at 1.
    Numpy stores coefficients in [x**2 + x + 1] order, so the fixed
    zeroth order denominator coefficent must comes last. (Edited.)
    """
    return np.polyval(p, x) / np.polyval(q, x)

def rational2_3(x, p0, p1, q0, q1, q2):
    return rational(x, [p0, p1], [q0, q1, q2])

def rational2_2(x, p0, p1, q0, q1):
    return rational(x, [p0, p1], [q0, q1])

x = np.linspace(0, 100, 300)  
y = (688.9*x + 114.4) / (x*x + 1038.0*x + 474.1)
#y = rational(x, [688.9, 114.4], [1.0, 1038.0, 474.1])
ynoise = y * (1.0 + np.random.normal(scale=0.01, size=x.shape))
popt, pcov = curve_fit(rational2_3, x, ynoise)
print popt

#plt.plot(x, y, label='original')
#plt.plot(x, ynoise, '.', label='data')
#plt.plot(x, rational2_3(x, *popt), label='fit')
#plt.legend()
#plt.show()


########### AWBS diffusive asymptotic #########################################
###############################################################################

###############################################################################
########### Spitzer-Harm diffusive asymptotic #################################
## Raw data used from SH_PR1953 paper.
vSH = np.array([0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.32, 0.34, 0.36, 0.38, 0.40, 0.44, 0.48, 0.52, 0.56, 0.60, 0.64, 0.68, 0.72, 0.76, 0.80, 0.88, 0.96, 1.04, 1.12, 1.20, 1.28, 1.36, 1.44, 1.52, 1.60, 1.76, 1.92, 2.08, 2.24, 2.40, 2.56, 2.72, 2.88, 3.04, 3.20])
dEZ1 = np.array([0.0008093, 0.001300, 0.001970, 0.002847, 0.003955, 0.005317, 0.006955, 0.008886, 0.01113, 0.01370, 0.01660, 0.02347, 0.0318, 0.04165, 0.05304, 0.06601, 0.08057, 0.09672, 0.1145, 0.1338, 0.1548, 0.2015, 0.2545, 0.3137, 0.3792, 0.4508, 0.5285, 0.6123, 0.7023, 0.7983, 0.9005, 1.123, 1.371, 1.645, 1.945, 2.273, 2.630, 3.017, 3.435, 3.887, 4.375, 5.465, 6.728, 8.190, 9.880, 11.83, 14.06, 16.62, 19.53, 22.74, 26.00])
dEZ100 = np.array([0.0001, 0.0001464, 0.0002074, 0.0002856, 0.0003842, 0.0005062, 0.0006554, 0.0008352, 0.001050, 0.001303, 0.001600, 0.002343, 0.003318, 0.004570, 0.006147, 0.008100, 0.01049, 0.01336, 0.01680, 0.02085, 0.02560, 0.03748, 0.05308, 0.07312, 0.09834, 0.1296, 0.1678, 0.2138, 0.2687, 0.3336, 0.4096, 0.5997, 0.8493, 1.170, 1.574, 2.074,2.684, 3.421, 4.300, 5.338, 6.554, 9.595, 13.59, 18.72, 25.18, 33.18, 42.95, 54.74, 68.80, 85.41, 104.9])
dTZ1 = np.array([0.0005906, 0.0009023, 0.001309, 0.001821, 0.002448, 0.003197, 0.004074, 0.005082, 0.006225, 0.007504, 0.008922, 0.01217, 0.01596, 0.02027, 0.02508, 0.03035, 0.03607, 0.04218, 0.04865, 0.05545, 0.06254, 0.07746, 0.09302, 0.1090, 0.1250, 0.1407, 0.1558, 0.1700, 0.1829, 0.1944, 0.2036, 0.2151, 0.2139, 0.1965, 0.1587, 0.0957, 0.0021, -0.1289, -0.3041, -0.5339, -0.8268, -1.657, -2.921, -4.774, -7.448, -11.23, -16.5, -23.5, -33.2, -46.0, -62.0])
dTZ100 = np.array([0.0001245, 0.0001821, 0.0002577, 0.0003546, 0.0004764, 0.0006271, 0.0008108, 0.001032, 0.001295, 0.001605, 0.001968, 0.002872, 0.004052, 0.005558, 0.007442, 0.009760, 0.01257, 0.01593, 0.01991, 0.02456, 0.02995, 0.04322, 0.06024, 0.08151, 0.1075, 0.1387, 0.1754, 0.2178, 0.2663, 0.3207, 0.3809, 0.5174, 0.6703, 0.8297, 0.9800, 1.099, 1.156, 1.113, 0.9167, 0.5060, -0.1966, -2.867, -8.061, -17.09, -31.69, -54.08, -87.05, -134.1, -199.3, -287.9, -405.8])
## Explicit values of Sh correction for Z = 1, 2, 4, 16, 10000
## obtained from corr = deltaT * epsilon / 0.4 in TABLE III
SH_Z = np.array([1, 2, 4, 16, 10000])
SH_gammaTovergammaE = np.array([0.46887895460797796, 0.6054441680081956, 0.7279908268569244, 0.8974525745257452, 1.0])
SH_corrs = np.array([0.23584070000000001, 0.36520749999999996, 0.5141982749999999, 0.7825953249999998, 1.0])
## Scale the velocity from v2Th to vth.
vSH *= 2.0**0.5
## Find a spline for the f0 data.
smooth = 0 # lower less smoothing
dEZ1_tck = splrep(vSH, dEZ1, s=smooth)
dEZ100_tck = splrep(vSH, dEZ100, s=smooth)
dTZ1_tck = splrep(vSH, dTZ1, s=smooth)
dTZ100_tck = splrep(vSH, dTZ100, s=smooth)
### SH Functions.
def SH_corr(Z):
    if (Z < SH_Z[0]):
        return SH_corrs[0]
    if (Z >= SH_Z[len(SH_Z)-1]):
        return SH_corrs[len(SH_Z)-1]
    for i in range(len(SH_Z)-1):
        if (SH_Z[i] <= Z < SH_Z[i+1]):
            return SH_corrs[i] + (SH_corrs[i+1] - SH_corrs[i]) / (SH_Z[i+1] - SH_Z[i]) * (Z - SH_Z[i])
def gammaTovergammaE(Z):
    if (Z < SH_Z[0]):
        return SH_gammaTovergammaE[0]
    if (Z >= SH_Z[len(SH_Z)-1]):
        return SH_gammaTovergammaE[len(SH_Z)-1]
    for i in range(len(SH_Z)-1):
        if (SH_Z[i] <= Z < SH_Z[i+1]):
            return SH_gammaTovergammaE[i] + (SH_gammaTovergammaE[i+1] - SH_gammaTovergammaE[i]) / (SH_Z[i+1] - SH_Z[i]) * (Z - SH_Z[i])
#print gammaTovergammaE(1.0), gammaTovergammaE(2.0), gammaTovergammaE(4.0), gammaTovergammaE(16.0), gammaTovergammaE(10000.0)
def SH_distribution(v, ne, Te, gradTe, Z, Gamma_ee):
    N = len(v)
    f1 = np.zeros(N) 
    j1 = np.zeros(N)
    q1 = np.zeros(N)
    ## Apply very rough interpolation between Z=1 and Z=100.
    dEZ1 = splev(v / vTh(Te), dEZ1_tck, der=0)
    dTZ1 = splev(v / vTh(Te), dTZ1_tck, der=0)
    dEZ100 = splev(v / vTh(Te), dEZ100_tck, der=0)
    dTZ100 = splev(v / vTh(Te), dTZ100_tck, der=0)
    dE = dEZ1 + (dEZ100 - dEZ1) / 99 * (Z - 1.0)
    dT = dTZ1 + (dTZ100 - dTZ1) / 99 * (Z - 1.0)
    D = 2.0 * dT + 3.0 / 2.0 * gammaTovergammaE(Z) * dE
    for i in range(N):
        vp = v[i]
        f1[i] = vTh(Te)**4.0 * 4.0 / Gamma_ee / ne / Z * D[i] * fM(ne, Te, vp) * gradTe / Te
        j1[i] = f1[i] * vp**3.0
        q1[i] = f1[i] * vp**5.0
    return f1, j1, q1
#dfdv = splev(vs, ftck, der=1)
########### Spitzer-Harm diffusive asymptotic #################################
###############################################################################

## Test of the oscillator equation.
if (0):
    N = 100
    #N = 10000
    v = np.linspace(0.0, 10.0, N)
    ## Prepare coefficients.
    a = np.ones(N)
    b = np.ones(N)
    c = np.ones(N)
    d = np.ones(N)
    a = 0.0*a*1e1
    b = b*1e0
    c = c*1e1
    d = d*0.0
    ## ICs.
    f0 = 1.0
    df0 = 0.0
    ## Solve second order ODE.
    f, df = ad2f_bdf_cf_equals_d(v, a, b, c, d, f0, df0)
    ff = bdf_cf_equals_d(v, b, c, d, f0)
    # Plot the result.
    plt.plot(v, f, 'b')
    plt.plot(v, ff, 'gx')
    plt.show()

def DistributionsOfZbar(N, ne, Te, dTedz, Z, G_ee):
    ## Number of cells in velocity magnitude.
    #N = 10000
    ## Multiple of thermal velocity for min(v).
    min_x_vTh = 0.001
    ## Multiple of thermal velocity for max(v).
    max_x_vTh = 7.0
    ## Define velocity magnitude discretization.
    vs = np.linspace(min_x_vTh * vTh(Te), max_x_vTh * vTh(Te), N)
    ## Define f1, j1, and q1 distributions of BGK. 
    fBGK1s = np.zeros(N)
    jBGK1s = np.zeros(N)
    qBGK1s = np.zeros(N)
    ## Lorentz gas E field.
    Ez_SH = vTh(Te)**2.0 * (dnedz / ne +  2.5 * dTedz / Te)
    ## Compute AWBS f1 part of the distribution.
    #corr = 2.0 # Heat flux fit
    #corr = 1.8 # Current fit
    ## Compute AWBS (you can check the zero current condition by varying Ez)
    #fAWBS1s, jAWBS1s, qAWBS1s = AWBS_distribution(vs, 0.0, ne, Te, dTedz, Z, Ez_SH, G_ee, corr)
    #fAWBS1s, jAWBS1s, qAWBS1s = _AWBS_distribution(vs, 0.0, ne, Te, dTedz, Z, Ez_SH, G_ee, corr)
    ## Generic values.
    corr_E = 1.0
    corr_nu = 0.5
    ## Exact values for Z = 1
    #corr_nu = 0.45642791622
    #corr_E = 0.999946330759
    fAWBS1s, jAWBS1s, qAWBS1s = AWBS_distribution(vs, ne, Te, dTedz, Z, Ez_SH, G_ee, corr_nu, corr_E)
    #corr = 1.7
    #fAWBS1s, jAWBS1s, qAWBS1s = HighVelocity_distribution(vs, 0.0, ne, Te, dTedz, Z, Ez_SH, G_ee, corr)
    ## Compute SH original distributions from the 1953 paper.
    fSH19531s, jSH19531s, qSH19531s = SH_distribution(vs, ne, Te, dTedz, Z, G_ee)
    ## Fill BGK distributions, which correspond to the Lorentz gas model, when
    ## scaling (Z + 0.24) / (Z + 4.2) * (Z + 1) / Z is applied
    corrSH = (Z + 0.24) / (Z + 4.2)
    print "(Z + 0.24) / (Z + 4.2): ", corrSH
    #corrSH = SH_corr(Z)
    for i in range(0, N):
        mfpei = vs[i]**4.0 / G_ee / ne / Z
        V2 = 0.5*vs[i]**2.0 / vTh(Te)**2.0
        fBGK1s[i] = - corrSH * mfpei * (V2 - 4.0) * fM(ne, Te, vs[i]) * dTedz / Te
        #fKIPP1s[i] = - mfpei * (3.0 / 8.0 * V2 - 3.0 / 4.0 / V2 - 1.0) * fM(ne, Te, vs[i]) * dTedz / Te
        jBGK1s[i] = fBGK1s[i] * vs[i]**3.0
        qBGK1s[i] = fBGK1s[i] * vs[i]**5.0
    return vs, fBGK1s, fAWBS1s, fSH19531s, jBGK1s, jAWBS1s, jSH19531s, qBGK1s, qAWBS1s, qSH19531s

## Perform actual distribution evaluation for given Zbar
N = 10000
## PRECISE EVALUATION OF Z dependence
FindAWBSdependenceOnZ(N, ne, Te, dTedz, Gamma_ee)
## PRECISE EVALUATION OF Z dependence
vs, fBGK1s, fAWBS1s, fSH19531s, jBGK1s, jAWBS1s, jSH19531s, qBGK1s, qAWBS1s, qSH19531s = DistributionsOfZbar(N, ne, Te, dTedz, Zbar, Gamma_ee)
## Perform rough evaluation for ZbarLv for BGK model.
vsLv_rough, fLvBGK1s_rough, fLvAWBS1s, fLvSH19531s, jLvBGK1s_rough, jLvAWBS1s, jLvSH19531s, qLvBGK1s_rough, qLvAWBS1s, qLvSH19531s = DistributionsOfZbar(30, ne, Te, dTedz, ZbarLv, GammaLv_ee)
## Perform detailed evaluation for ZbarLv for AWBS model. 
vsLv, fLvBGK1s, fLvAWBS1s, fLvSH19531s, jLvBGK1s, jLvAWBS1s, jLvSH19531s, qLvBGK1s, qLvAWBS1s, qLvSH19531s = DistributionsOfZbar(N, ne, Te, dTedz, ZbarLv, GammaLv_ee)

## Integrate currents and fluxes.
jBGK = sum(jBGK1s)
#jKIPP = sum(jKIPP1s)
jAWBS = sum(jAWBS1s)
jSH1953 = sum(jSH19531s)
print "jBGK, jAWBS, jSH1953: ", jBGK, jAWBS, jSH1953
qBGK = sum(qBGK1s)
#qKIPP = sum(qKIPP1s)
qAWBS = sum(qAWBS1s)
qSH1953 = sum(qSH19531s)
print "qBGK, qAWBS, qSH1953: ", qBGK, qAWBS, qSH1953
print "|qBGK - qAWBS|/ qBGK: ", abs(qBGK - qAWBS) / qBGK 
qLvAWBS = sum(qLvAWBS1s)
qLvBGK = sum(qLvBGK1s)
print "|qLvBGK - qLvAWBS|/ qLvBGK: ", abs(qLvBGK - qLvAWBS) / qLvBGK
## Plot SH and KIPP normalized distribution functions.
vs_norm = vs / vTh(Te)
#vs_norm = vs / (2.0**0.5*vTh(Te))
#fig, ax1 = plt.subplots()
#ax1.plot(vs_norm, fBGK1s, 'r', label=r'$f_1-$SH')
#ax1.plot(vs_norm, fKIPP1s, 'b', label=r'$f_1-$KIPP')
#ax1.legend(loc='upper right', fancybox=True, framealpha=0.8)
#ax1.set_title('Zbar '+str(Zbar))
#plt.show()
fig, ax1 = plt.subplots()
ax1.plot(vs_norm, jBGK1s, 'g-.', label=r'$j_1-$BGK')
#ax1.plot(vs_norm, jKIPP1s, 'b', label=r'$j_1-$KIPP')
ax1.plot(vs_norm, jAWBS1s, 'r--', label=r'$j_1-$AWBS')
ax1.plot(vs_norm, jSH19531s, 'k', label=r'$j_1-$FP')
ax1.legend(loc='upper right', fancybox=True, framealpha=0.8)
ax1.set_title('Zbar '+str(Zbar))
plt.show()
fig, ax1 = plt.subplots()
ax1.plot(vs_norm, qBGK1s, 'g-.', label=r'$q_1^{Z=1}$BGK')
ax1.plot(vs_norm, qAWBS1s, 'r--', label=r'$q_1^{Z=1}$AWBS')
ax1.plot(vs_norm, qSH19531s, 'k', label=r'$q_1^{Z=1}$SH')
ax1.legend(loc='upper left', fancybox=True, framealpha=0.8)
ax1.set_xlabel(r'$v / v_{th}$')
ax1.set_ylabel(r'$q_1^{Z=1}$ [a.u.]')
ax2 = ax1.twinx()
ax2.plot(vs_norm, qLvAWBS1s, 'y', label=r'$q_1^{Z=116}$AWBS')
ax2.plot(vsLv_rough / vTh(Te), qLvBGK1s_rough, 'bx', label=r'$q_1^{Z=116}$Lorentz')
ax2.legend(loc='upper right', fancybox=True, framealpha=0.8)
ax2.set_ylabel(r'$q_1^{Z=116}$ [a.u.]')
ax1.set_title('Distribution flux moment')
fig.tight_layout()
for ext in ["png", "pdf", "eps"]:
   print("saving q1s.%s" % (ext,))
   plt.savefig("q1s.%s" % (ext,), bbox_inches="tight")
plt.show()

###############################################################################
########### AWBS kinetic numerical experiments results ########################

## marker size.
ms = 8

Kn_AladinZ2 = np.array([2e-4, 3.5e-3, 4.1e-2, 9.2e-2, 1.7e-1, 2.8e-1, 5.6e-1])
q_AladinZ2 = np.array([3.25e13, 1.99e13, 1.023e13, 4.35e12, 2.79e12, 2.13e12, 10.4e11])
qSH_AladinZ2 = np.array([3.26e13, 2.2e13, 1.828e13, 1.13e13, 9.82e12, 8.98e12, 7.15e12])
qC7_AladinZ2 = np.array([3.2e13, 1.97e13, 9.11e12, 3.22e12, 1.69e12, 8.46e11, 1.32e11])

plt.plot(np.log10(Kn_AladinZ2), qC7_AladinZ2 / qSH_AladinZ2, 'b*', markersize=ms, label=r'$q^{Z=2}$AP1')
plt.plot(np.log10(Kn_AladinZ2), q_AladinZ2 / qSH_AladinZ2, 'rv', markersize=ms, label=r'$q^{Z=2}$Aladin')
#plt.plot(np.log10(Kn_AladinZ2), np.log10(qC7_AladinZ2 / qSH_AladinZ2), 'b*', markersize=ms, label=r'$q^{Z=2}$AP1')
#plt.plot(np.log10(Kn_AladinZ2), np.log10(q_AladinZ2 / qSH_AladinZ2), 'rv', markersize=ms, label=r'$q^{Z=2}$Aladin')

Kn_Impact = np.array([4.9e-2, 3.8e-2, 1.9e-2, 1.3e-3])
q_Impact = np.array([1.104e13, 9.8e12, 6.45e12, 6.98e11])
qSH_Impact = np.array([2.11e13, 1.72e13, 8.93e12, 7.05e11])
qC7_Impact = np.array([9.78e12, 8.82e12, 6.25e12, 6.88e11])

plt.plot(np.log10(Kn_Impact), q_Impact / qSH_Impact, 'go', markersize=ms, label=r'$q^{Z=2}$Impact')
plt.plot(np.log10(Kn_Impact), qC7_Impact / qSH_Impact, 'b*', markersize=ms)
#plt.plot(np.log10(Kn_Impact), np.log10(q_Impact / qSH_Impact), 'go', markersize=ms, label=r'$q^{Z=2}$Impact')
#plt.plot(np.log10(Kn_Impact), np.log10(qC7_Impact / qSH_Impact), 'b*', markersize=ms)

Kn_Calder = np.array([3.5e-2])
q_Calder = np.array([1.06e13])
qSH_Calder = np.array([1.77e13])
qC7_Calder = np.array([9.2e12])

plt.plot(np.log10(Kn_Calder), q_Calder / qSH_Calder, 'ks', markersize=ms, label=r'$q^{Z=2}$Calder')
plt.plot(np.log10(Kn_Calder), qC7_Calder / qSH_Calder, 'b*', markersize=ms)
#plt.plot(np.log10(Kn_Calder), np.log10(q_Calder / qSH_Calder), 'ks', markersize=ms, label=r'$q^{Z=2}$Calder')
#plt.plot(np.log10(Kn_Calder), np.log10(qC7_Calder / qSH_Calder), 'b*', markersize=ms)

Kn_AladinZ10 = np.array([5.1e-3, 3.5e-2])
q_AladinZ10 = np.array([5.63e12, 1.36e12])
qSH_AladinZ10 = np.array([9.23e12, 5.06e12])
qC7_AladinZ10 = np.array([5.19e12, 6.97e11])

#plt.plot(np.log10(Kn_AladinZ10), q_AladinZ10 / qSH_AladinZ10, 'rx', markersize=ms, label=r'$q^{Z=10}$Aladin')
#plt.plot(np.log10(Kn_AladinZ10), qC7_AladinZ10 / qSH_AladinZ10, 'bx', markersize=ms, label=r'$q^{Z=10}$C7')

plt.xlabel(r'$\log_{10}$(Kn$^e$)')
plt.ylabel(r'$q/q_{SH}$')
#plt.ylabel(r'$\log_{10}(q/q_{SH})$')
plt.legend(loc='lower left', fancybox=True, framealpha=0.8)
plt.tight_layout()
for ext in ["png", "pdf", "eps"]:
   print("saving Kn_results.%s" % (ext,))
   plt.savefig("Kn_results.%s" % (ext,), bbox_inches="tight")
plt.show()

########### AWBS kinetic numerical experiments results ########################
###############################################################################

quit()

## Check the validity of action of FP collision operator on fSH.
ad2fdv2_SH, bdfdv_SH, cf_SH, d_SH = AllFokkerPlanckEquationTerms(vs, fMs, fBGK1s, dfMdzs, Ez_SH, Gamma_ee)
## Check the validity of action of FP collision operator on fKIPP.
ad2fdv2_KIPP, bdfdv_KIPP, cf_KIPP, d_KIPP = AllFokkerPlanckEquationTerms(vs, fMs, fKIPP1s, dfMdzs, Ez_SH, Gamma_ee)
## Plot the terms of FP action.
vs_norm = vs / vTh(Te)
## Distribution fSH.
fig, ax1 = plt.subplots()
ax1.plot(vs_norm, ad2fdv2_SH, 'r', label='ad2fdv2')
ax1.plot(vs_norm, bdfdv_SH, 'b', label='bdfdv')
ax1.plot(vs_norm, cf_SH, 'g', label='cf')
ax1.plot(vs_norm, d_SH, 'k', label='d')
ax1.legend(loc='upper right', fancybox=True, framealpha=0.8)
ax1.set_title('fSH - Zbar '+str(Zbar))
plt.show()
## Distribution fKIPP.
fig, ax1 = plt.subplots()
ax1.plot(vs_norm, ad2fdv2_KIPP, 'r', label='ad2fdv2')
ax1.plot(vs_norm, bdfdv_KIPP, 'b', label='bdfdv')
ax1.plot(vs_norm, cf_KIPP, 'g', label='cf')
ax1.plot(vs_norm, d_KIPP, 'k', label='d')
ax1.legend(loc='upper right', fancybox=True, framealpha=0.8)
ax1.set_title('fKIPP - Zbar '+str(Zbar))
plt.show()

quit()

## In order to start the iteration of f1, set fSH1 as initial state.
f1s = fSH19531s

## The diffusive electron Fokker-Planck equation coefficients for f1.
## These remain fixed during the iteration/convergence of f1.
a, b, ce, ci = FPcoefficientsF0(vs, fMs)
c = 0.0 * ce + ci

vs_norm = vs / vTh(Te)
fig, ax1 = plt.subplots()
ax1.plot(vs_norm, a, 'r', label='a')
ax1.plot(vs_norm, b, 'b', label='b')
ax1.plot(vs_norm, c, 'g', label='c')
ax1.legend(loc='upper right', fancybox=True, framealpha=0.8)
ax1.set_title('Zbar '+str(Zbar))
plt.show()

## Right hand side of the diffusive electron Fokker-Planck equation. 
d2f0dv2If1_df0dvIf1, df0dz, Edf0dv = FPtermsF0F1(vs, fMs, f1s, dfMdzs, Ez_SH, Gamma_ee)
d = d2f0dv2If1_df0dvIf1 + df0dz + Edf0dv

## Analytic RHS computed directly from fM function properties. Sanity check.
d2fMdv2If1_dfMdvIf1, dfMdz, EdfMdv = FPtermsFMF1(vs, fMs, f1s, dfMdzs, Ez_SH, Te, Gamma_ee)
dM = 0.0 * d2fMdv2If1_dfMdvIf1 + dfMdz + EdfMdv

## All the terms of the diffusive electron Fokker-Planck equation in order to
## analyze the importance of the derivatives of f1 with respect to v, 
## which are usually omitted in FP codes.
f = f1s
f = 0.0 * f
for i in range(5):
   ad2fdv2, bdfdv, cf, d = AllFokkerPlanckEquationTerms(vs, fMs, f, dfMdzs, Ez_SH, Gamma_ee)

   a = a * 0.0
   b = b * 0.0
   #d = d - ad2fdv2 - bdfdv
   f = bdf_cf_equals_d(vs, b, c, dM, 0.0)
   #f, df = ad2f_bdf_cf_equals_d(vs, a, b, c, d, 0.0, 0.0)

   vs_norm = vs / vTh(Te)
   fig, ax1 = plt.subplots()
   ax1.set_title('Zbar '+str(Zbar))
   ax1.plot(vs_norm, fMs, 'r', label='fM')
   ax1.legend(loc='upper left')
   ax2 = ax1.twinx()
   ax2.plot(vs_norm, fBGK1s, 'b', label='fBGK1')
   ax2.plot(vs_norm, f, 'b--', label='f1')
   ax2.legend(loc='upper right')
   plt.show()

## See the scattering dominance of on ions with respect to Zbar.
vs_norm = vs / vTh(Te)
fig, ax1 = plt.subplots()
ax1.plot(vs_norm, cf, 'r', label='(ci + ce) * f1')
ax1.plot(vs_norm, ci * f1s, 'r:', label='ci * f1')
ax1.plot(vs_norm, bdfdv, 'b', label='C1(f0) * df1dv')
ax1.plot(vs_norm, ad2fdv2, 'b--', label='C2(f0) * d2f1dv2')
ax1.legend(loc='upper right', fancybox=True, framealpha=0.8)
ax1.set_title('Zbar '+str(Zbar))
plt.show()

vs_norm = vs / vTh(Te)
fig, ax1 = plt.subplots()
ax1.plot(vs_norm, d, 'r', label='rhs')
ax1.plot(vs_norm, df0dz + Edf0dv, 'r--', label='v*df0dz + Edf0dv')
ax1.legend(loc='upper right', fancybox=True, framealpha=0.8)
ax1.set_title('Zbar '+str(Zbar))
plt.show()

#### FP equation diffusive regime #############################################
###############################################################################

###############################################################################
########### AWBS diffusive asymptotic ######################################### 
def solve_bweuler(v, f0, ne, Te, gradTe, Z, E):
    N = len(v)
    f1 = np.zeros(N) 
    f1[0] = f0
    for i in range(N-1):
        dv = v[i+1] - v[i]
        vp = v[i]
        rhs = Z/vp*((vp**2.0/2.0/vTh(Te)**2.0 - 1.5)*gradTe/Te - E/vTh(Te)**2.0)
        rhs = rhs * fM(ne, Te, vp)
        f1[i+1] = (f1[i] + dv*rhs)/(1.0 + dv*(3.0 - Z)/vp) # ee isotropization
        #f1[i+1] = (f1[i] + dv*rhs)/(1.0 + dv*(4.0 - Z)/vp)
    return f1

## Multiples of thermal velocity setting the velocity space range.
ml_max = 10.0
ml_min = 0.05
## The heat flux after integration takes the form
## qH = me/Zbar/sigma/coulLog*128/(2*pi)**0.5*(kB/me)**(7/2)*T**(5/2)*gradT,
## where mfp_ei = v**4/sigma/coulLog/ni/Zbar/Zbar, 
## i.e. sigma corresponds to unit charges collisions.
corr = (688.9*Zbar + 114.4)/(Zbar**2.0 + 1038*Zbar + 474.1)
#print "Zbar, corr:", Zbar, corr
cmag = 1./corr
## Classical Lorentz approximation electric field.
xi = 2.5
## SH correction.
#xi = 1.0 + 1.5 * (Zbar + 0.477) / (Zbar + 2.15)
Efield = vTh(Te)**2.0 * (dnedz / ne +  xi * dTedz / Te)

## Solve the AWBS ODE problem.
N = 5000
v = np.linspace(ml_max*vTh(Te), ml_min*vTh(Te), N)
dv = (ml_max - ml_min)*vTh(Te)/(N-1)
## Use both explicit and implicit method of ODE solve.
#sol = solve_bweuler(v, 0.0, ne, Te, dTedz, Zbar, Efield)
cmag = 0.5
fBGK1s = solve_bweuler(v, 0.0, ne, Te, dTedz, cmag*Zbar, Efield)

"""
###############################################################################
########### Kinetic integration of AWBS and C7 ################################ 
###############################################################################
fM_analytic = np.zeros(N)
SHf1 = np.zeros(N)
SHj = np.zeros(N)
SHq = np.zeros(N)
SHJ = 0.0
SHQ = 0.0
AWBSf1_corr = np.zeros(N)
AWBSj_corr = np.zeros(N)
AWBSq_corr = np.zeros(N)
AWBSJ_corr = 0.0
AWBSQ_corr = 0.0
AWBSf1 = np.zeros(N)
AWBSj = np.zeros(N)
AWBSq = np.zeros(N)
AWBSJ = 0.0
AWBSQ = 0.0
for i in range(N):
    vp = v[i]
    ## The mean free path has standard v^4 dependence, sigma is cross section
    ## given by model and Zbar increases the effect of Coulomb potential in
    ## ei collisions.
    mfp_ei = vp**4.0/sigma/coulLog/ni/Zbar/Zbar
    SHf1[i] = - (Zbar + 0.24)/(Zbar + 4.2)*((vp**2.0/2.0/vTh(Te)**2.0 - 4.0)*gradTe/Te)*fM(vp, Te)*vp*vp
    #SHf1[i] = - (Zbar + 0.24)/(Zbar + 4.2)*((vp**2.0/2.0/vTh(Te)**2.0 - 1.5)*gradTe/Te - Efield/vTh(Te)**2.0)*fM(vp, Te)*vp*vp
    fM_analytic[i] = fM(vp, Te)*vp*vp
    SHj[i] = mfp_ei*vp*SHf1[i]
    SHq[i] = mfp_ei*me/2.0*vp*vp*vp*SHf1[i]
    SHJ = SHJ + 4.0*pi/3.0*SHj[i]*dv
    SHQ = SHQ + 4.0*pi/3.0*SHq[i]*dv
    AWBSf1_corr[i] = sol_corr[i]*vp*vp
    AWBSj_corr[i] = mfp_ei*vp*AWBSf1_corr[i]
    AWBSq_corr[i] = mfp_ei*vp*vp*vp*me/2.0*AWBSf1_corr[i]
    AWBSJ_corr = AWBSJ_corr + 4.0*pi/3.0*AWBSj_corr[i]*dv
    AWBSQ_corr = AWBSQ_corr + 4.0*pi/3.0*AWBSq_corr[i]*dv
    AWBSf1[i] = sol[i]*vp*vp
    AWBSj[i] = mfp_ei*vp*AWBSf1[i]
    AWBSq[i] = mfp_ei*vp*vp*vp*me/2.0*AWBSf1[i]
    AWBSJ = AWBSJ + 4.0*pi/3.0*AWBSj[i]*dv
    AWBSQ = AWBSQ + 4.0*pi/3.0*AWBSq[i]*dv

#######################################
## Analytic SH formula ################
#######################################
## Analytical formula from AWBShos.pdf, providing the Lorentz gas approximation
## further multiplied by SH low Z factor.
mfp_ei = (vTh(Te))**4.0/sigma/coulLog/ni/Zbar/Zbar
## We add one to Zbar (i.e. Zbar+1.) in order to include the ee collisions.
mfp_tot = (vTh(Te))**4.0/sigma/coulLog/ni/Zbar/(Zbar+1.)
## Classical length scale definition.
##L = Te / abs(gradTe)
## More accurate length scale definition.
L = dTe / abs(gradTe)
## SH nu correction.
SHcorr = (Zbar + 0.24)/(Zbar + 4.2)
## SH Efield correction.
SHcorr = SHcorr * (3.5 - xi)
## SH analytic formula of heat flux.
SHQ_analytic = - SHcorr * 128.0/(2.0*pi)**0.5*ne*vTh(Te)*kB*Te*mfp_ei*gradTe/Te
## Standar formulation used in hydrocodes.
SHQ_pete = - SHcorr * 1.31e10 / coulLog / Zbar * Te**2.5 * gradTe
## Apply variation of sigma.
SHQ_pete = SHQ_pete * 8.1027575e17 / sigma
Kn =  mfp_tot / L
Kn_flux = SHQ_analytic / (SHcorr * 128.0/(2.0*pi)**0.5 * ne * vTh(Te) * kB * Te)
Kn_pete = SHQ_pete / (SHcorr * 128.0/(2.0*pi)**0.5 * ne * vTh(Te) * kB * Te)
## Express flux proportionality with respect to SHQ_analytic.
proporC7EQ = C7EQ / SHQ_analytic

#######################################
## Print comparison results ###########
#######################################
## Show the Knudsen number
print 'Kn: ', Kn, 'mfp_tot[microns]: ', mfp_tot*1e4
print 'Kn_flux, Kn_pete: ', Kn_flux, Kn_pete
## Print integrated values
print "SHQ[erg/s/cm2]:              ", SHQ
print "SHQ_analytic[erg/s/cm2]:     ", SHQ_analytic
print "SHQ_pete[erg/s/cm2]:         ", SHQ_pete
print "C7EQ[erg/s/cm2]:             ", C7EQ
print "diffAWBSQ_corr[erg/s/cm2]:   ", AWBSQ_corr
print "diffAWBSQ[erg/s/cm2]:        ", AWBSQ
#print "SHJ:        ", SHJ
#print "AWBSJ_corr: ", AWBSJ_corr
#print "AWBSJ:      ", AWBSJ

###############################################################################
########### Plotting of results ###############################################
###############################################################################
###############################################################################
import matplotlib.pyplot as plt
import matplotlib
## Global setting of plotting.
font = {'family' : 'Sans',
        #'weight' : 'bold',
        'size'   : 18}
figure = {'figsize' : '10.5, 6.5'} # a tuple in inches
matplotlib.rc('font', **font)
matplotlib.rc('figure', **figure)

###############################################################################
########### Kinetic SH, Diffusive AWBS and C7 comparison ######################
###############################################################################
## Shrink the x axis appropriately.
mult = 8
p_v = v[v < mult*vTh(Te)]
p_fM_analytic = fM_analytic[v < mult*vTh(Te)]
p_SHq = SHq[v < mult*vTh(Te)]
p_AWBSq_corr = AWBSq_corr[v < mult*vTh(Te)]
p_AWBSq = AWBSq[v < mult*vTh(Te)]
if (args.C7):
   p_C7Ev = C7Ev[C7Ev < mult*vTh(Te)]
   p_C7Emehalff1v5 = C7Emehalff1v5[C7Ev < mult*vTh(Te)]
   p_C7Emehalff0v5 = C7Emehalff0v5[C7Ev < mult*vTh(Te)]
   p_C7Ef0v2 = C7Emehalff0v5[C7Ev < mult*vTh(Te)]
   ## out-of-equilibrium source
   for i in range(len(p_C7Ev)):
      p_C7Ef0v2[i] = fM(p_C7Ev[i], Te)*p_C7Ev[i]*p_C7Ev[i] + p_C7Emehalff0v5[i] / (4.0*pi) / me * 2.0 / p_C7Ev[i] / p_C7Ev[i] / p_C7Ev[i]
## Set labels.
fig, ax1 = plt.subplots()
ax1.set_ylabel(r'$q_1 = m_e v^2/2\, v f_1 v^2$ [a.u.]')
ax1.set_xlabel('v/vT')
print "ne: ", ne
print "Kn: ", Kn
ax1.set_title('Kinetics (Z='+"{:.1f}".format(float(Zbar))+r', n$_e$='+"{:.1e}".format(float(ne))+', Kn='+"{:.1e}".format(Kn)+')')
## Plot kinetic analysis.
if (args.kinSH):
   ax1.plot(p_v/vTh(Te), p_SHq, SHcolor+"-", label=r'$q_1^{SH}$')
if (AWBSoriginal):
   ax1.plot(p_v/vTh(Te), p_AWBSq, "b-.", label=r'$q_1^{AWBS}$')
if (AWBSstar):
   ax1.plot(p_v/vTh(Te), p_SHq * (3.0/8.0*p_v*p_v/vTh(Te)/vTh(Te) - 3.0*vTh(Te)*vTh(Te)/p_v/p_v - 2.0)/(p_v*p_v/vTh(Te)/vTh(Te)-8.0), "g--", label=r'$q_1^{KIPP}$')
   #ax1.plot(p_v/vTh(Te), p_AWBSq_corr, "r"+"-.", label=r'$q_1^{AWBS^*}$')
if (args.C7):
   ax1.plot(p_C7Ev/vTh(Te), p_C7Emehalff1v5 / (4.0*pi/3.0), C7Ecolor+'-', label=r'$q_1^{C7}$')
   #ax1.plot(p_C7Ev/vTh(Te), p_C7Emehalff1v5 / (4.0*pi/3.0), C7Ecolor+'-', label=r'$q_1^{C7E}$'+'('+"{:.2f}".format(proporC7EQ)+r'$q_h^{SH}$)')
## q0 axis
ax2 = ax1.twinx()
if (args.kinSH):
   ax2.plot(p_v/vTh(Te), p_fM_analytic, SHcolor+':', label=r'$f_0^{SH}$')
   #ax2.plot(p_v/vTh(Te), me / 2.0 * p_v * p_v * p_v * p_fM_analytic, SHcolor+':', label=r'$q_0^{SH}$')
if (args.C7):
   ax2.plot(p_C7Ev/vTh(Te), p_C7Ef0v2, C7Ecolor+'--', label=r'$f_0^{C7}$')
   #ax2.plot(p_C7Ev/vTh(Te), p_C7Emehalff0v5 / (4.0*pi) / me * 2.0 / p_C7Ev / p_C7Ev / p_C7Ev, C7Ecolor+'--', label=r'$f_0^{C7}$')
   #ax2.plot(p_C7Ev/vTh(Te), p_C7Emehalff0v5 / (4.0*pi), C7Ecolor+'--', label=r'$q_0^{C7E}$')
ax2.set_ylabel(r'$f_0 v^2$ [a.u.]')
#ax2.set_ylabel(r'$q_0 = m_e v^2/2\, v f_0 v^2$ [a.u.]')
ax1.legend(loc='upper right', fancybox=True, framealpha=0.8)
ax2.legend(loc='lower right', fancybox=True, framealpha=0.8)
for ext in ["png", "pdf", "eps"]:
   print("saving kinetics.%s" % (ext,))
   plt.savefig("kinetics.%s" % (ext,), bbox_inches="tight")
if (args.pltshow):
   plt.show()

if (0):
   fig, ax1 = plt.subplots()
   ax1.set_ylabel(r'$j_1 =  q_e v f_1 v^2$ [a.u.]')
   ax1.set_xlabel('v/vT')
   ax1.set_title('Kinetics (Z='+str(Zbar)+r', n$_e$='+"{:.1e}".format(ne)+', Kn='+"{:.1e}".format(Kn)+')')
   ## Plot kinetic analysis.
   ax1.plot(p_v/vTh(Te), p_SHq / p_v / p_v / me * qe, SHcolor+"-", label=r'$j_1^{SH}$')
   if (AWBSoriginal):
      ax1.plot(p_v/vTh(Te), p_AWBSq / p_v / p_v / me * qe, "b-.", label=r'$j_1^{AWBS}$')
   if (AWBSstar):
      ax1.plot(p_v/vTh(Te), p_SHq * (3.0/8.0*p_v*p_v/vTh(Te)/vTh(Te) - 3.0*vTh(Te)*vTh(Te)/p_v/p_v - 2.0)/(p_v*p_v/vTh(Te)/vTh(Te)-8.0) / p_v / p_v / me * qe, "g--", label=r'$j_1^{KIPP}$')
   if (args.C7):
      ax1.plot(p_C7Ev/vTh(Te), p_C7Emehalff1v5 / (4.0*pi/3.0) / p_C7Ev / p_C7Ev / me * qe, C7Ecolor+'-', label=r'$j_1^{C7}$')
   ax1.legend(loc='upper right', fancybox=True, framealpha=0.8)
   for ext in ["png", "pdf", "eps"]:
      print("saving j_kinetics.%s" % (ext,))
      plt.savefig("j_kinetics.%s" % (ext,), bbox_inches="tight")
   if (args.pltshow):
      plt.show()
"""

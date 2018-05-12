import numpy as np
from math import pi
from math import exp
import matplotlib.pyplot as plt
import matplotlib

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


## Some default values.
ne = 1.0e23
Te = 100.0
dTedz = -1e0
dnedz = 0.0
Zbar = 1000.0

## Full electron-electron cross section.
Gamma_ee = sigma * coulLog
#Gamma_ee = 1.0

lambda_ei = vTh(Te)**4.0 / Gamma_ee / ne / (Zbar + 1.0)
Kn_ei = lambda_ei * dTedz / Te
print "Kn_ei: ", Kn_ei 

import argparse
## Create parser object.
parser = argparse.ArgumentParser(description='Compare the diffusive asymptotic to C7 computation.')
## Define input arguments.
parser.add_argument("-s", "--sigma", help="Sigma for electro-ion cross-section.", type=float)
parser.add_argument("-Np", "--Nproc", help="Number of processors used to compute the data.", type=int)
## A no value argument solution.
parser.add_argument("-ps", "--pltshow", action='store_true', help="Plot show() by adding -ps/--pltshow argument.")
parser.add_argument("-lF1", "--labelFluxExt1", help="Use -lF1/--labelFluxExt1 to use and label VFPdata/flux1.dat.")


## Parse arguments.
args = parser.parse_args()
#if (args.sigma):
#   sigma = args.sigma

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

def ad2f_bdf_cf_equals_d(vs, a, b, c, d, f0, df0):
    ## Solving equation a*dfdv2 + b*dfdv + c*f = d
    ## omega = sqrt(a / c), theta = b / 2 / sqrt(a * c)
    ## Create discretization of f.
    omega = (a / c)**0.5
    theta = b / 2.0 / (a * c)**0.5
    ## Find number of cells and use equidistant velocity step.
    N = len(vs)
    dv =  (max(vs) - min(vs)) / (N - 1)
    print "omega, theta, omega/dv: ", omega[0], theta[0], omega[0] / dv
    ## Define the unknowns.
    f = np.zeros(N)
    df = np.zeros(N)
    ### Set initial conditions for 0->infinity integration.
    #f[0] = f0
    #df[0] = df0
    ## Perform backward Euler integration.
    #for i in range(N-1):
    #    aa = (a[i] + a[i+1]) / 2.0
    #    bb = (b[i] + b[i+1]) / 2.0
    #    cc = (c[i] + c[i+1]) / 2.0
    #    dd = (d[i] + d[i+1]) / 2.0
    #    df[i+1] = (dv * dd + aa * df[i] - dv * cc * f[i]) / (dv**2.0 * cc + aa + dv * bb)
    #    f[i+1] = (dv * dd + aa * df[i] - (aa + dv * bb) * df[i+1]) / dv / cc
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

## Test of the oscillator equation.
if (1):
    N = 10000
    v = np.linspace(0.0, 10.0, N)
    ## Prepare coefficients.
    a = np.ones(N)
    b = np.ones(N)
    c = np.ones(N)
    d = np.ones(N)
    a = a*1e1
    b = b*1e0
    c = c*1e1
    d = d*0.0
    ## ICs.
    f0 = 1.0
    df0 = 0.0
    ## Solve second order ODE.
    f, df = ad2f_bdf_cf_equals_d(v, a, b, c, d, f0, df0)
    # Plot the result.
    plt.plot(v, f)
    plt.show()

## Number of cells in velocity magnitude.
N = 2000
## Multiple of thermal velocity for min(v).
min_x_vTh = 0.1
## Multiple of thermal velocity for max(v).
max_x_vTh = 7.0
## Define velocity magnitude discretization.
vs = np.linspace(min_x_vTh * vTh(Te), max_x_vTh * vTh(Te), N)
## Define f0 and f1 discretization.
fMs = np.zeros(N)
dfMdzs = np.zeros(N)
f1s = np.zeros(N)
## Fill the Maxwell-Boltzmann distribution discretizations.
for i in range(0, N):
    fMs[i] = fM(ne, Te, vs[i])
    dfMdzs[i] = dfMdz(ne, Te, dnedz, dTedz, vs[i])
## Fill f1 with Lorentz gas approximation as a starting value.
for i in range(0, N):
    mfpei = vs[i]**4.0 / Gamma_ee / ne / (Zbar + 1.0)
    f1s[i] = - mfpei * (0.5*vs[i]**2.0 / vTh(Te)**2.0 - 4.0) * fM(ne, Te, vs[i]) * dTedz / Te
## Lorentz gas E field.
Ez = vTh(Te)**2.0 * (dnedz / ne +  2.5 * dTedz / Te)

## Integrate the f0-Rosenbluth potentials. 
## These remain constant during the f1 iteration/convergence.
I0, I2, Jm1 = RosenbluthPotentialsF0(vs, fMs)

I1, I3, Jm2 = RosenbluthPotentialsF1(vs, f1s)

#print "fMs"
#print fMs
#print "f1s"
#print f1s

#print "I0"
#print I0
#print "I2"
#print I2
#print "Jm1"
#print Jm1

#print "I1"
#print I1
#print "I3"
#print I3
#print "Jm2"
#print Jm2

a_fM = (I2 + Jm1) / 3.0 / vs

b_fM = (3.0 * I0 - I2 + 2.0 * Jm1) / 3.0 / vs**2.0

ce_fM = 8.0*pi*fMs - (3.0 * I0 - I2 + 2.0 * Jm1) / 3.0 / vs**3.0

ci = - Zbar * ne / vs**3.0

d_fM = 1.0 / Gamma_ee * (vs * dfMdzs - vs * Ez / vTh(Te)**2.0 * fMs)

d_f1 = - fMs * 1.0 / 15.0 / vs / vTh(Te)**2.0 * (3.0 * vs**2.0 / vTh(Te)**2.0 * (I3 + Jm2) - 5.0 * (I1 + Jm2))

#print "a_fM"
#print a_fM
#print "b_fM"
#print b_fM
#print "ce_fM"
#print ce_fM
#print "ci"
#print ci
#print "d_fM"
#print d_fM
#print "d_f1"
#print d_f1
#print "d_f10"
#print d_f10

#from scipy.interpolate import splev, splrep
#smooth = 0 # lower less smoothing
### Find a spline for the f1 Lorentz data.
#f1Ltck = splrep(vs, f1s, s=smooth)
#df1Ldv = splev(vs, f1Ltck, der=1)
#d2f1Ldv2 = splev(vs, f1Ltck, der=2)
### Find a spline for the fM data.
#fMtck = splrep(vs, fMs, s=smooth)
#dfMdv = splev(vs, fMtck, der=1)
#d2fMdv2 = splev(vs, fMtck, der=2)
#
# Numerical check of d_f1.
#d_f10 = - (1.0 / 5.0 / vs * (I3 + Jm2) * d2fMdv2 + 1.0 / 15.0 / vs**2.0 * (5.0 * I1 - 3.0 * I3 + 2.0 * Jm2) * dfMdv)
#
#print "a_fM * d2f1Ldv2"
#print a_fM * d2f1Ldv2
#print "b_fM * df1Ldv"
#print b_fM * df1Ldv
#print "ce_fM * f1"
#print ce_fM * f1s
#print "ci * f1"
#print ci * f1s

## Global setting of plotting.
font = {'family' : 'Sans',
        #'weight' : 'bold',
        'size'   : 18}
figure = {'figsize' : '10.5, 6.5'} # a tuple in inches
matplotlib.rc('font', **font)
matplotlib.rc('figure', **figure)

plt.plot(vs, fMs, label='fM')
plt.plot(vs, f1s, label='f1')
plt.legend()
plt.title('Zbar '+str(Zbar))
plt.show()

plt.plot(vs, d_f1, label='d_f1')
plt.plot(vs, d_fM, label='d_fM')
plt.legend()
plt.title('Zbar '+str(Zbar))
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
sol = solve_bweuler(v, 0.0, ne, Te, dTedz, Zbar, Efield)
sol_corr = solve_bweuler(v, 0.0, ne, Te, dTedz, cmag*Zbar, Efield)

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

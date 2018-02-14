import numpy as np
from math import pi
from math import exp

## Fundamental physical constants in cgs. 
kB = 1.3807e-16
me = 9.1094e-28

## Number of processors used to run C7.
Nproc = 8
maxProcOrder = 6 # Corresponds to maximum of million processors.

## Some default values.
ne = 5.0e19
Te = 10000.0
gradTe = -1.0
Zbar = 4.0
sigma = 1e15
xpoint = 0.5
Ecorrect=False
Emimic=False

import argparse
## Create parser object.
parser = argparse.ArgumentParser(description='Analyze the diffusive asymptotic of AWBS and compare to C7 computation.')
## Define input arguments.
parser.add_argument("-n", "--ne", help="Electron density at the point.", type=float)
parser.add_argument("-s", "--sigma", help="Electro-ion cross-section.", type=float)
parser.add_argument("-Z", "--Zbar", help="Ionization at the point.", type=float)
parser.add_argument("-xp", "--xpoint", help="Kinetic analysis at this point.", type=float)
parser.add_argument("-Np", "--Nproc", help="Number of processors used to compute the data.", type=int)
## A no value argument solution.
parser.add_argument("-Er", "--Ecorrect", action='store_true', help="Display the Ecorrect computation results.")
parser.add_argument("-Em", "--Emimic", action='store_true', help="Display the Emimic computation results.")

## Parse arguments.
args = parser.parse_args()
if args.ne:
    ne = args.ne
if args.sigma:
    sigma = args.sigma
if args.Zbar:
    Zbar = args.Zbar
if args.xpoint:
    xpoint = args.xpoint
if args.Nproc:
    Nproc = args.Nproc
if args.Ecorrect:
    Ecorrect = args.Ecorrect
if args.Emimic:
    Emimic = args.Emimic
###############################################################################
########### Loading of results of parallel C7 code ############################
###############################################################################
###############################################################################
## Load the profiles provided by C7.
def loadC7data(Nproc, file_base):
   ## Global lists.
   C7x_raw = []
   C7rho_raw = []
   C7Te_raw = []
   C7intf0_raw = []
   C7j_raw = []
   C7q_raw = []
   ## Gather together the lists from each processors output.
   print "Loading data from files:"
   for proc in range(Nproc):
      file = file_base+str(proc).zfill(maxProcOrder)
      C7x_proc, C7rho_proc, C7Te_proc, C7intf0_proc, C7j_proc, C7q_proc = np.loadtxt(file,  usecols=(0, 1, 2, 3, 4, 5), unpack=True)
      C7x_raw.extend(C7x_proc)
      C7rho_raw.extend(C7rho_proc)
      C7Te_raw.extend(C7Te_proc)
      C7intf0_raw.extend(C7intf0_proc)
      C7j_raw.extend(C7j_proc)
      C7q_raw.extend(C7q_proc)
      print file
   ## Sort the lists with respect to the position x.
   sorted_indices = np.array(C7x_raw).argsort()
   C7x = np.array([C7x_raw[sorted_indices[j]] for j in range(len(C7x_raw))])
   C7rho = np.array([C7rho_raw[sorted_indices[j]] for j in range(len(C7x_raw))])
   C7Te = np.array([C7Te_raw[sorted_indices[j]] for j in range(len(C7x_raw))])
   C7intf0 = np.array([C7intf0_raw[sorted_indices[j]] for j in range(len(C7x_raw))])
   C7j = np.array([C7j_raw[sorted_indices[j]] for j in range(len(C7x_raw))])
   C7q = np.array([C7q_raw[sorted_indices[j]] for j in range(len(C7x_raw))])
   return C7x, C7rho, C7Te, C7intf0, C7j, C7q

C7x, C7rho, C7Te, C7intf0_Ec, C7j_Ec, C7q_Ec = loadC7data(Nproc, 'Ecorrect_data/C7_1_profiles.')
C7x, C7rho, C7Te, C7intf0_Em, C7j_Em, C7q_Em = loadC7data(Nproc, 'Emimic_data/C7_1_profiles.')
###############################################################################
########### Analysis of diffusive asymptotic of AWBS model #################### 
###############################################################################
###############################################################################
from scipy.interpolate import splev, splrep 
## Obtain a given point value and gradient of temperature.
smooth = 0 # lower less smoothing
## Find a spline for the temperature data.
tck = splrep(C7x, C7Te, s=smooth)
## Evaluate at the given point xpoint.
xp = np.array(xpoint)
## Assign temperature profile values for further analysis.
Te = splev(xp, tck, der=0)
gradTe = splev(xp, tck, der=1)
print "xpoint, ne, sigma, Zbar, Te, gradTe: ", xpoint, ne, sigma, Zbar, Te, gradTe
## Assign a whole profile of SH flux.
C7gradTe = splev(C7x, tck, der=1)

## Ion density is used as reference.
ni = ne / Zbar

###############################################################################
########### AWBS diffusive asymptotic ######################################### 
###############################################################################
def vTh(T): 
    return (kB*T/me)**0.5
def fM(v, T):
    return ne/(vTh(T)**3.0*(2.0*pi)**1.5)*exp(-v**2.0/2.0/vTh(T)**2.0)
def solve_bweuler(v, f0, T, gradT, Z, E):
    N = len(v)
    f1 = np.zeros(N) 
    f1[0] = f0
    for i in range(N-1):
        dv = v[i+1] - v[i]
        vp = v[i]
        rhs = Z/vp*((vp**2.0/2.0/vTh(T)**2.0 - 1.5)*gradT/T - E/vTh(T)**2.0)
        rhs = rhs * fM(vp, T)
        f1[i+1] = (f1[i] + dv*rhs)/(1.0 + dv*(4.0 - Z)/vp)
    return f1

## Multiples of thermal velocity setting the velocity space range.
ml_max = 10.0
ml_min = 0.05
## The heat flux after integration takes the form
## qH = me/Zbar/sigma*128/(2*pi)**0.5*(kB/me)**(7/2)*T**(5/2)*gradT,
## where mfp_ei = v**4/sigma/ni/Zbar, i.e. sigma corresponds to ee collisions.
corr = (688.9*Zbar + 114.4)/(Zbar**2.0 + 1038*Zbar + 474.1)
#print "Zbar, corr:", Zbar, corr
cmag = 1./corr
## Classical Lorentz approximation electric field.
Efield = vTh(Te)**2.0*2.5*gradTe/Te

## Solve the AWBS ODE problem.
N = 5000
v = np.linspace(ml_max*vTh(Te), ml_min*vTh(Te), N)
dv = (ml_max - ml_min)*vTh(Te)/(N-1)
## Use both explicit and implicit method of ODE solve.
sol = solve_bweuler(v, 0.0, Te, gradTe, Zbar, Efield)
sol_corr = solve_bweuler(v, 0.0, Te, gradTe, cmag*Zbar, Efield)

###############################################################################
########### Kinetic integration of AWBS and C7 ################################ 
###############################################################################
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
    mfp_ei = vp**4.0/sigma/ni/Zbar
    SHf1[i] = - (Zbar + 0.24)/(Zbar + 4.2)*((vp**2.0/2.0/vTh(Te)**2.0 - 1.5)*gradTe/Te - Efield/vTh(Te)**2.0)*fM(vp, Te)*vp*vp
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
## Load C7 kinetic results ############
#######################################
## No explicit treatment of Efield, we use mimicing by reducing ne in source.
C7v, C7mehalff1v5 = np.loadtxt('Emimic_data/fe_point_Emimic.txt',  usecols=(0, 4), unpack=True)
C7Q = 0.0
NC7 = C7v.size - 1
for i in range(NC7):
    dC7v = C7v[i] - C7v[i+1]
    C7Q = C7Q + C7mehalff1v5[i]*dC7v
## Explicit treatment of Efield.
C7Ev, C7Emehalff1v5 = np.loadtxt('Ecorrect_data/fe_point_Ecorrect.txt',  usecols=(0, 4), unpack=True)
C7EQ = 0.0
NC7E = C7Ev.size - 1
for i in range(NC7E):
    dC7Ev = C7Ev[i] - C7Ev[i+1]
    C7EQ = C7EQ + C7Emehalff1v5[i]*dC7Ev

#######################################
## Analytic SH formula ################
#######################################
## Analytical formula from AWBShos.pdf, providing the Lorentz gas approximation
## further multiplied by SH low Z factor.
mfp_ei = (vTh(Te))**4.0/sigma/ni/Zbar
L = Te / abs(gradTe)
SHQ_analytic = - (Zbar + 0.24)/(Zbar + 4.2) * 128.0/(2.0*pi)**0.5*ne*vTh(Te)*kB*Te*mfp_ei*gradTe/Te
Kn = mfp_ei/L
Kn_flux = SHQ_analytic / ((Zbar + 0.24)/(Zbar + 4.2) * 128.0/(2.0*pi)**0.5 * ne * vTh(Te) * kB * Te)
## Express flux proportionality with respect to SHQ_analytic.
proporC7EQ = C7EQ / SHQ_analytic
proporC7Q = C7Q / SHQ_analytic

#######################################
## C7 SH flux profile #################
#######################################
## Calculate the whole profile of SH flux according to Te from C7.
C7SHQ_analytic = - (Zbar + 0.24)/(Zbar + 4.2) * 128.0/(2.0*pi)**0.5*ne*vTh(C7Te)*kB*C7Te*(vTh(C7Te))**4.0/sigma/ni/Zbar*C7gradTe/C7Te

#######################################
## Print comparison results ###########
#######################################
## Show the Knudsen number
print 'Kn: ', Kn, 'Kn from flux: ', Kn_flux 
## Print integrated values
print "SHQ:              ", SHQ
print "SHQ_analytic:     ", SHQ_analytic
print "C7EQ:             ", C7EQ
print "C7Q:              ", C7Q
print "diffAWBSQ_corr:   ", AWBSQ_corr
print "diffAWBSQ:        ", AWBSQ
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

## Given line styles.
lsC7 = 'k--'
lsC7E = 'b-.'
lsAWBS = 'r--'
lsSH = 'g:'
lblC7 = r'C7$^*$'
lblC7E = r'C7E'
lblAWBS = r'AWBS$^*$'
lblSH = r'SH'

###############################################################################
########### C7 plasma profiles ################################################
###############################################################################
# It is useful to plot the profiles with respect to microns.
C7x_microns = np.array(C7x) * 1e4
## Set labels.
#plt.xlabel(r'z [$\mu$m]')
#plt.ylabel(r'$\rho$ [g/cm$^3$]')
#plt.title(r'Density')
#plt.plot(C7x_microns, C7rho)
#plt.show()
## Set labels.
plt.xlabel(r'z [$\mu$m]')
plt.ylabel(r'$T_e$ [eV]')
plt.title(r'Electron temperature')
plt.plot(C7x_microns, C7Te)
plt.show()
## Set labels.
#plt.xlabel(r'z [$\mu$m]')
#plt.ylabel(r'$\int f_0 dv$ [1/cm$^3$]')
#plt.title(r'Isotropic part of f')
#if (Ecorrect):
#   plt.plot(C7x_microns, C7intf0_Ec, lsC7E, label=lblC7E)
#if (Emimic):
#   plt.plot(C7x_microns, C7intf0_Em, lsC7, label=lblC7)
#plt.show()
## Set labels.
plt.xlabel(r'z [$\mu$m]')
plt.ylabel(r'$j$ [1/s/cm$^2$]')
plt.title(r'Flux/current')
if (Ecorrect):
   plt.plot(C7x_microns, C7j_Ec, lsC7E, label=lblC7E)
if (Emimic):
   plt.plot(C7x_microns, C7j_Em, lsC7, label=lblC7)
plt.legend()
plt.show()
## Set labels.
plt.xlabel(r'z [$\mu$m]')
plt.ylabel(r'q [W/cm$^2$]')
plt.title(r'Heat flux (Kn='+"{:.1e}".format(Kn)+')')
plt.plot(C7x_microns, C7SHQ_analytic, lsSH, label=lblSH)
if (Ecorrect):
   plt.plot(C7x_microns, C7q_Ec, lsC7E, label=lblC7E)
if (Emimic):
   plt.plot(C7x_microns, C7q_Em, lsC7, label=lblC7)
plt.legend()
plt.show()

###############################################################################
########### Kinetic SH, Diffusive AWBS and C7 comparison ######################
###############################################################################
## Shrink the x axis appropriately.
mult = 8
p_v = v[v < mult*vTh(Te)]
p_SHq = SHq[v < mult*vTh(Te)]
p_AWBSq_corr = AWBSq_corr[v < mult*vTh(Te)]
p_AWBSq = AWBSq[v < mult*vTh(Te)]
p_C7v = C7v[C7v < mult*vTh(Te)]
p_C7mehalff1v5 = C7mehalff1v5[C7v < mult*vTh(Te)]
p_C7Ev = C7Ev[C7Ev < mult*vTh(Te)]
p_C7Emehalff1v5 = C7Emehalff1v5[C7Ev < mult*vTh(Te)]


## Set labels.
plt.ylabel(r'$q_1 = m_e v^2/2\, v f_1 v^2$ [a.u.]')
plt.xlabel('v/vT')
plt.title('Z = '+str(Zbar)+', Kn = '+"{:.1e}".format(Kn))
## Plot kinetic analysis.
plt.plot(p_v/vTh(Te), p_SHq, lsSH, label=lblSH)
#plt.plot(p_v/vTh(Te), p_AWBSq, 'g-.', label='AWBS')
plt.plot(p_v/vTh(Te), p_AWBSq_corr, 'r--', label=lblAWBS)
if (Ecorrect):
   if (len(C7Ev)<30):
      plt.plot(p_C7Ev/vTh(Te), p_C7Emehalff1v5 / (4.0*pi/3.0), 'bx', label=lblC7E+'('+"{:.2f}".format(proporC7EQ)+r'q$_{SH}$)')
   else:
      plt.plot(p_C7Ev/vTh(Te), p_C7Emehalff1v5 / (4.0*pi/3.0), lsC7E, label=lblC7E+'('+"{:.2f}".format(proporC7EQ)+r'q$_{SH}$)')
#plt.plot(C7v/vTh(Te), 1.5 * C7mehalff1v5 / (4.0*pi/3.0), 'k:', label=r'C7')
if (Emimic):
   plt.plot(p_C7v/vTh(Te), p_C7mehalff1v5 / (4.0*pi/3.0), lsC7, label=lblC7+'('+"{:.2f}".format(proporC7Q)+r'q$_{SH}$)')
plt.legend(loc='best')
plt.show()

"""
###############################################################################
########### Obtain and plot Z-dependent AWBS correction #######################
###############################################################################
## Physical fix by a magic constant.
## Z = 1
#cmag = 1.882
## Z = 5
#cmag = 1.59575
## Z = 10
#cmag = 1.565
## Z = 20
#cmag = 1.56
## Z = 50
#cmag = 1.585
## Z = 100
#cmag = 1.65
## Z = 200
#cmag = 1.8
coeffs = np.array([0.5313496280552604, 0.6266645777847407, 0.6389776357827476, 0.641025641025641, 0.6309148264984227, 0.6060606060606061, 0.5555555555555556])
Zbars = [1.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0]
p1 = 688.9
p2 = 114.4
q1 = 1038
q2 = 474.1
N = 1000
ZZ = np.linspace(min(Zbars), max(Zbars), N)
icoeffs = np.array([(p1*ZZ[i] + p2)/(ZZ[i]**2.0 + q1*ZZ[i] + q2) for i in range(N)])
##print icoeffs
#Zbar = Zbars[0]
#corr = coeffs[0]
#Zbar = Zbars[1]
#corr = coeffs[1]
#Zbar = Zbars[2]
#corr = coeffs[2]
#Zbar = Zbars[3]
#corr = coeffs[3]
#Zbar = Zbars[4]
#corr = coeffs[4]
#Zbar = Zbars[5]
#corr = coeffs[5]
#Zbar = Zbars[6]
#corr = coeffs[6]

#plt.plot(Zbars, coeffs, 'rx', label='pointwise corrections')
#plt.plot(ZZ, icoeffs, 'b', label=str(p1)+r'$Z + $'+str(p2)+r'$/Z^2 + $'+str(q1)+r'$Z + $'+str(q2))
#plt.legend(loc='best')
#plt.xlabel('Z')
#plt.title('Rational function fit')
"""

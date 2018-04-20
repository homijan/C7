import numpy as np
from math import pi
from math import exp

## Fundamental physical constants in cgs. 
kB = 1.6022e-12
me = 9.1094e-28
qe = 4.8032e-10

## Number of processors used to run C7.
Nproc = 8
maxProcOrder = 6 # Corresponds to maximum of million processors.

## Coulomb logarithm.
coulLog = 10

## Some default values.
ni = 1.0e20
Te = 1000.0
gradTe = -1.0
Zbar = 4.0
sigma = 8.1027575e17 ## Matching the SH diffusive flux.
AWBSstar=False
AWBSoriginal=False
usexpoint=False
vlimshow=False

import argparse
## Create parser object.
parser = argparse.ArgumentParser(description='Compare the diffusive asymptotic to C7 computation.')
## Define input arguments.
parser.add_argument("-cl", "--coulLog", help="Coulomb logarithm for electro-ion cross-section.", type=float)
parser.add_argument("-Np", "--Nproc", help="Number of processors used to compute the data.", type=int)
## A no value argument solution.
parser.add_argument("-ps", "--pltshow", action='store_true', help="Plot show() by adding -ps/--pltshow argument.")
parser.add_argument("-As", "--AWBSstar", action='store_true', help="Display the AWBS* diffusive asymptotic by adding -As/--AWBSstar argument.")
parser.add_argument("-Ao", "--AWBSoriginal", action='store_true', help="Display the AWBSoriginal diffusive asymptotic by adding -Ao/--AWBSoriginal argument.")
parser.add_argument("-C7", "--C7", action='store_true', help="Display the C7 computation results by adding -C7/--C7 argument.")
parser.add_argument("-xp", "--usexpoint", action='store_true', help="Use an xpoint computed output of distribution function instead of nonlocal heat flux maximum.")
parser.add_argument("-vs", "--vlimshow", action='store_true', help="vlim show by adding -vs/--vlimshow argument.")
## Faking arguments providing different labels than C7E and C7*.
#parser.add_argument("-lC7", "--labelC7", help="Force to use -lC7/--labelC7 label for C7_data.")

## Parse arguments.
args = parser.parse_args()
if (args.coulLog):
   coulLog = args.coulLog

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
   C7j_raw = []
   C7Ex_raw = []
   C7q_raw = []
   C7corrE_raw = []
   C7ne_raw = []
   C7zbar_raw = []
   ## Gather together the lists from each processors output.
   print "Loading data from files:"
   for proc in range(Nproc):
      file = file_base+str(proc).zfill(maxProcOrder)
      C7x_proc, C7rho_proc, C7Te_proc, C7j_proc, C7Ex_proc, C7q_proc, C7corrE_proc, C7ne_proc, C7zbar_proc = np.loadtxt(file,  usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8), unpack=True)
      C7x_raw.extend(C7x_proc)
      C7rho_raw.extend(C7rho_proc)
      C7Te_raw.extend(C7Te_proc)
      C7j_raw.extend(C7j_proc)
      C7Ex_raw.extend(C7Ex_proc)
      C7q_raw.extend(C7q_proc)
      C7corrE_raw.extend(C7corrE_proc)
      C7ne_raw.extend(C7ne_proc)
      C7zbar_raw.extend(C7zbar_proc)
      print file
   ## Sort the lists with respect to the position x.
   sorted_indices = np.array(C7x_raw).argsort()
   C7x = np.array([C7x_raw[sorted_indices[j]] for j in range(len(C7x_raw))])
   C7rho = np.array([C7rho_raw[sorted_indices[j]] for j in range(len(C7x_raw))])
   C7Te = np.array([C7Te_raw[sorted_indices[j]] for j in range(len(C7x_raw))])
   C7j = np.array([C7j_raw[sorted_indices[j]] for j in range(len(C7x_raw))])
   C7Ex = np.array([C7Ex_raw[sorted_indices[j]] for j in range(len(C7x_raw))])
   C7q = np.array([C7q_raw[sorted_indices[j]] for j in range(len(C7x_raw))])
   C7corrE = np.array([C7corrE_raw[sorted_indices[j]] for j in range(len(C7x_raw))])
   C7ne = np.array([C7ne_raw[sorted_indices[j]] for j in range(len(C7x_raw))])
   C7zbar = np.array([C7zbar_raw[sorted_indices[j]] for j in range(len(C7x_raw))])
   return C7x, C7rho, C7Te, C7j, C7Ex, C7q, C7corrE, C7ne, C7zbar

## Load C7 data.
C7x, C7rho, C7Te, C7j, C7Ex, C7q, C7corrE, C7ne, C7zbar = loadC7data(Nproc, 'C7_data/C7_1_profiles.')
###############################################################################
########### Analysis of diffusive asymptotic of AWBS model #################### 
###############################################################################
###############################################################################
from scipy.interpolate import splev, splrep 
## Obtain a given point value and gradient of temperature.
smooth = 0 # lower less smoothing
## Find a spline for the temperature data.
Tetck = splrep(C7x, C7Te, s=smooth)
## Assign a whole temperature profile.
C7gradTe = splev(C7x, Tetck, der=1)
## Find a spline for the electron density data.
netck = splrep(C7x, C7ne, s=smooth)
## Assign a whole electron density profile.
C7gradne = splev(C7x, netck, der=1)
## Find a spline for the ionization data.
zbartck = splrep(C7x, C7zbar, s=smooth)
## Assign a whole ionization profile.
C7zbar = splev(C7x, zbartck, der=0)

#######################################
## Load C7 kinetic results ############
#######################################
## Explicit treatment of Efield.
if (usexpoint):
   C7Expoints, C7Ev, C7Emehalff1v5, C7Emehalff0v5 = np.loadtxt('C7_data/fe_point_C7.txt',  usecols=(0, 1, 5, 6), unpack=True)
   C7Expoint = C7Expoints[0]
   C7EQ = 0.0
   NC7E = C7Ev.size - 1
   for i in range(NC7E):
      dC7Ev = C7Ev[i] - C7Ev[i+1] 
      C7EQ = C7EQ + C7Emehalff1v5[i]*dC7Ev
else:
   C7Expoints, C7Ev, C7Emehalff1v5, C7Emehalff0v5 = np.loadtxt('C7_data/fe_pointmax_C7.txt',  usecols=(0, 1, 5, 6), unpack=True)
   C7Expoint = C7Expoints[0]
   C7EQ = 0.0
   NC7E = C7Ev.size - 1
   for i in range(NC7E):
      dC7Ev = C7Ev[i] - C7Ev[i+1]
      C7EQ = C7EQ + C7Emehalff1v5[i]*dC7Ev

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
        f1[i+1] = (f1[i] + dv*rhs)/(1.0 + dv*(3.0 - Z)/vp) # ee isotropization
        #f1[i+1] = (f1[i] + dv*rhs)/(1.0 + dv*(4.0 - Z)/vp)
    return f1

## Use the C7E xpoint as reference.
xpoint = C7Expoint
## Evaluate at the given point xpoint.
xp = np.array(xpoint)
## Assign temperature, electron density and ionization profile values 
## for a further analysis.
Te = splev(xp, Tetck, der=0)
gradTe = splev(xp, Tetck, der=1)
ne = splev(xp, netck, der=0)
Zbar = splev(xp, zbartck, der=0)
## Ion density.
ni = ne / Zbar
## For more accurate Kn evaluate the "differential" of temperature.
mfp_ei = vTh(Te)**4.0/sigma/coulLog/ni/Zbar/Zbar
Te0 = C7Te.min()
dTe = abs(Te - Te0)
print "xpoint, ni, ne, sigma, coulLog, Zbar, Te, gradTe, dTe: ", xpoint, ni, ne, sigma, coulLog, Zbar, Te, gradTe, dTe

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
Efield = vTh(Te)**2.0*xi*gradTe/Te

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
Kn =  mfp_tot / L
Kn_flux = SHQ_analytic / (SHcorr * 128.0/(2.0*pi)**0.5 * ne * vTh(Te) * kB * Te)
Kn_pete = SHQ_pete / (SHcorr * 128.0/(2.0*pi)**0.5 * ne * vTh(Te) * kB * Te)
## Express flux proportionality with respect to SHQ_analytic.
proporC7EQ = C7EQ / SHQ_analytic


#######################################
## C7 SH flux profile #################
#######################################
## Calculate the whole profile of SH flux according to Te from C7.
## We add one to Zbar (i.e. Zbar+1.) in order to include the ee collisions.
C7mfp_ei = (vTh(C7Te))**4.0/sigma/coulLog/ni/Zbar/(Zbar+1.)
Temin_reference = 0.9*C7Te.min()
C7Kn = C7mfp_ei * abs(C7gradTe) / (C7Te - Temin_reference)
## SH nu correction.
C7SHcorr = (C7zbar + 0.24)/(C7zbar + 4.2)
C7SHQ_analytic = - C7SHcorr * 128.0/(2.0*pi)**0.5*ne*vTh(C7Te)*kB*C7Te*(vTh(C7Te))**4.0/sigma/coulLog/ni/C7zbar/C7zbar*C7gradTe/C7Te
C7SHE_analytic = vTh(C7Te)**2.0*(C7gradne/C7ne + xi*C7gradTe/C7Te)

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

## Given line styles.
lsC7 = 'k--'
lsC7E = 'b-.'
lsAWBSc = 'r--'
lsAWBSo = 'r:'
lsSH = 'g:'
lblC7 = r'C7*'
lblC7E = r'C7E'
lblAWBSc = r'AWBS*'
lblAWBSo = r'AWBS'
lblSH = r'SH'
## Force to use an explicit labels, if passed.
#if args.labelC7:
#    lblC7E = args.labelC7

###############################################################################
########### C7 plasma profiles ################################################
###############################################################################
# It is useful to plot the profiles with respect to microns.
C7x_microns = np.array(C7x) * 1e4

## r=2 corresponds better
#xmicronsSH, QWcm2SH = np.loadtxt('../../VFPdata/flux_SH.dat',  usecols=(0, 1), unpack=True)
#xmicronsVFP, QWcm2VFP = np.loadtxt('../../VFPdata/flux_VFP.dat',  usecols=(0, 1), unpack=True)
#xmicronsSNB, QWcm2SNB = np.loadtxt('../../VFPdata/flux_SNB.dat',  usecols=(0, 1), unpack=True)

SHcolor = 'k'
C7Ecolor = 'r'
## Set labels.
fig, ax1 = plt.subplots()
ax1.set_xlabel(r'z [$\mu$m]')
ax1.set_ylabel(r'$q_h$ [W/cm$^2$]'+r', $T_e\in$('+"{:.0f}".format(C7Te.min())+', '+"{:.0f}".format(C7Te.max())+') [eV]')
ax1.set_title(r'Heat flux (Z = '+"{:.1f}".format(float(Zbar))+r', $\lambda_{th}$='+"{:.4f}".format(mfp_tot*1e4)+r'[$\mu$m])')
#ax1.set_title(r'Heat flux (Z = '+str(Zbar)+', Kn='+"{:.1e}".format(Kn)+')')
## Heat fluxes are displayed in W/cm2, i.e. energy is converted from ergs to J.
#ax1.plot(xmicronsSH, QWcm2SH, 'y-', label=r'SH')
#ax1.plot(xmicronsVFP, QWcm2VFP, 'y-', label=r'$q_h^{VFP}$')
#ax1.plot(xmicronsSNB, QWcm2SNB, 'm-', label=r'$q_h^{SNB}$')

ax1.plot(C7x_microns, C7SHQ_analytic * 1e-7, SHcolor+'-', label=r'$q_h^{SH}$')
if (args.C7):
   ax1.plot(C7x_microns, C7q * 1e-7, C7Ecolor+'-', label=r'$q_h^{C7E}$')
## Special treatment of temperature profile.
C7Te_scaled = C7Te*(C7SHQ_analytic.max() - C7SHQ_analytic.min())/(C7Te.max() - C7Te.min())
C7Te_scaled = C7Te_scaled - (C7Te_scaled.max() - C7SHQ_analytic.max())
ax1.plot(C7x_microns, C7Te_scaled * 1e-7, 'b:', label=r'$T_e$')
## Second Efield axis.
ax2 = ax1.twinx()
if (vlimshow):
   ax2.set_ylabel(r'E [a.u.]'+r', $v_{lim}/v_{th}\in$('+"{:.1f}".format(C7corrE.min())+', '+"{:.0f}".format(C7corrE.max())+')' )
else:
   ax2.set_ylabel(r'E [a.u.]')
ax2.plot(C7x_microns, me/qe*C7SHE_analytic, SHcolor+'-.', label=r'E$^{SH}$')
if (args.C7):
   ax2.plot(C7x_microns, me/qe*C7Ex, C7Ecolor+'--', label=r'E$^{C7E}$')
## Special treatment of the corrE showing the limit velocity/vTh 
## to be affected by E field.
if (vlimshow):
   C7corrE_scaled = C7corrE*(C7Ex.max() - C7Ex.min())/(C7corrE.max() - C7corrE.min()) 
   C7corrE_scaled = C7corrE_scaled - (C7corrE_scaled.max() - C7Ex.max())
   ax2.plot(C7x_microns, me/qe*C7corrE_scaled, 'k-.', label=r'$v_{lim}/v_{th}$')
fig.tight_layout()
ax1.legend(loc='center left', fancybox=True, framealpha=0.8)
ax2.legend(loc='center right', fancybox=True, framealpha=0.8)
for ext in ["png", "pdf", "eps"]:
   print("saving heatflux.%s" % (ext,))
   plt.savefig("heatflux.%s" % (ext,), bbox_inches="tight")
if (args.pltshow):
   plt.show()

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

## Set labels.
fig, ax1 = plt.subplots()
ax1.set_ylabel(r'$q_1 = m_e v^2/2\, v f_1 v^2$ [a.u.]')
ax1.set_xlabel('v/vT')
print "ne: ", ne
print "Kn: ", Kn
ax1.set_title('Kinetics (Z='+"{:.1f}".format(float(Zbar))+r', n$_e$='+"{:.1e}".format(float(ne))+', Kn='+"{:.1e}".format(Kn)+')')
## Plot kinetic analysis.
ax1.plot(p_v/vTh(Te), p_SHq, SHcolor+"-", label=r'$q_1^{SH}$')
if (AWBSoriginal):
   ax1.plot(p_v/vTh(Te), p_AWBSq, "b-.", label=r'$q_1^{AWBS}$')
if (AWBSstar):
   ax1.plot(p_v/vTh(Te), p_SHq * (3.0/8.0*p_v*p_v/vTh(Te)/vTh(Te) - 3.0*vTh(Te)*vTh(Te)/p_v/p_v - 2.0)/(p_v*p_v/vTh(Te)/vTh(Te)-8.0), "g--", label=r'$q_1^{KIPP}$')
   #ax1.plot(p_v/vTh(Te), p_AWBSq_corr, "r"+"-.", label=r'$q_1^{AWBS^*}$')
if (args.C7):
   ax1.plot(p_C7Ev/vTh(Te), p_C7Emehalff1v5 / (4.0*pi/3.0), C7Ecolor+'-', label=r'$q_1^{C7E}$')
   #ax1.plot(p_C7Ev/vTh(Te), p_C7Emehalff1v5 / (4.0*pi/3.0), C7Ecolor+'-', label=r'$q_1^{C7E}$'+'('+"{:.2f}".format(proporC7EQ)+r'$q_h^{SH}$)')
## q0 axis
ax2 = ax1.twinx()
ax2.plot(p_v/vTh(Te), p_fM_analytic, SHcolor+':', label=r'$f_0^{SH}$')
#ax2.plot(p_v/vTh(Te), me / 2.0 * p_v * p_v * p_v * p_fM_analytic, SHcolor+':', label=r'$q_0^{SH}$')
if (args.C7):
   ax2.plot(p_C7Ev/vTh(Te), p_C7Emehalff0v5 / (4.0*pi) / me * 2.0 / p_C7Ev / p_C7Ev / p_C7Ev, C7Ecolor+'--', label=r'$f_0^{C7E}$')
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
import numpy as np
from math import pi
from math import exp

## Physical constants in cgs.
kB = 1.3807e-16
me = 9.1094e-28
## Default values. 
ne = 1.0e20
sigma = 1e17
Zbar = 4.0

import argparse
## Create parser object.
parser = argparse.ArgumentParser(description='Analyze diffusive asymptotic of AWBS + compare to C7.')
## Define input arguments.
parser.add_argument("-n", "--ne", help="Electron density at the point.", type=float)
parser.add_argument("-s", "--sigma", help="Electro-ion cross-section.", type=float)
parser.add_argument("-Z", "--Zbar", help="Ionization at the point.", type=float)
parser.add_argument("-xp", "--xpoint", help="Kinetic analysis at this point.", type=float)

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

## Load the profiles provided by C7.
## Number of processors used to run C7.
Nproc = 8
maxProcOrder = 6 # Corresponds to maximum of million processors.
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
    file = '../tmp/C7_1_profiles.'+str(proc).zfill(maxProcOrder)
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

"""
# We use ion density as reference.
ni = ne / Zbar

def vTh(T): 
    return (kB*T/me)**0.5

def fM(v, T):
    return ne/(vTh(T)**3.0*(2.0*pi)**1.5)*exp(-v**2.0/2.0/vTh(T)**2.0)

ml_max = 10.0
ml_min = 0.05
# The heat flux after integration takes the form
# qH = me/Zbar/sigma*128/(2*pi)**0.5*(kB/me)**(7/2)*T**(5/2)*gradT,
# where mfp_ei = v**4/sigma/ni/Zbar, i.e. sigma corresponds to ei collisions.
Efield = vTh(Te)**2.0*2.5*gradTe/Te

N = 5000
v = np.linspace(ml_max*vTh(Te), ml_min*vTh(Te), N)
dv = (ml_max - ml_min)*vTh(Te)/(N-1)

# Post-process transport values
SHf1 = np.zeros(N)
SHj = np.zeros(N)
SHq = np.zeros(N)
SHJ = 0.0
SHQ = 0.0
for i in range(N):
    vp = v[i]
    # The mean free path has standard v^4 dependence, sigma is cross section
    # given by model and Zbar increases the effect of Coulomb potential in
    # ei collisions
    mfp_ei = vp**4.0/sigma/ni/Zbar
    SHf1[i] = - (Zbar + 0.24)/(Zbar + 4.2)*((vp**2.0/2.0/vTh(Te)**2.0 - 1.5)*gradTe/Te - Efield/vTh(Te)**2.0)*fM(vp, Te)*vp*vp
    SHj[i] = mfp_ei*vp*SHf1[i]
    SHq[i] = mfp_ei*me/2.0*vp*vp*vp*SHf1[i]
    SHJ = SHJ + 4.0*pi/3.0*SHj[i]*dv
    SHQ = SHQ + 4.0*pi/3.0*SHq[i]*dv

# Analytical formula from AWBShos.pdf, providing the Lorentz gas approximation
# further multiplied by SH low Z factor.
mfp_ei = (vTh(Te))**4.0/sigma/ni/Zbar
L = Te / abs(gradTe)
SHQ_analytic = - (Zbar + 0.24)/(Zbar + 4.2) * 128.0/(2.0*pi)**0.5*ne*vTh(Te)*kB*Te*mfp_ei*gradTe/Te
Kn = mfp_ei/L
Kn_flux = SHQ_analytic / ((Zbar + 0.24)/(Zbar + 4.2) * 128.0/(2.0*pi)**0.5 * ne * vTh(Te) * kB * Te)
# Show the Knudsen number
print 'Kn: ', Kn, 'Kn from flux: ', Kn_flux 

# Print integrated values
print "SHQ:          ", SHQ
print "SHQ_analytic: ", SHQ_analytic
"""

import matplotlib.pyplot as plt
import matplotlib
font = {'family' : 'Sans',
        #'weight' : 'bold',
        'size'   : 16}
matplotlib.rc('font', **font)

### Shrink the x axis appropriately.
#mult = 8
#p_v = v[v < mult*vTh(Te)]
#p_SHq = SHq[v < mult*vTh(Te)]
#
### Set labels.
#plt.ylabel(r'$q_1 = m_e/2 v^2 f_1 v^2$ [a.u.]')
#plt.xlabel('v/vT')
#plt.title('Z = '+str(Zbar)+', Kn = '+"{:.1e}".format(Kn))
### Plot kinetic analysis.
#plt.plot(p_v/vTh(Te), p_SHq, 'b', label='SH')
#plt.legend(loc='best')
#plt.grid()
#plt.show()

# It is useful to plot the profiles with respect to microns.
C7x_microns = C7x * 1e4
## Set labels.
plt.xlabel(r'z [$\mu$m]')
plt.ylabel(r'$\rho$ [g/cm$^3$]')
plt.title(r'Density')
plt.plot(C7x_microns, C7rho)
plt.show()
## Set labels.
plt.xlabel(r'z [$\mu$m]')
plt.ylabel(r'$T_e$ [eV]')
plt.title(r'Electron temperature')
plt.plot(C7x_microns, C7Te)
plt.show()
## Set labels.
plt.xlabel(r'z [$\mu$m]')
plt.ylabel(r'$\int f_0 dv$ [1/cm$^3$]')
plt.title(r'Isotropic part of f')
plt.plot(C7x_microns, C7intf0)
plt.show()
## Set labels.
plt.xlabel(r'z [$\mu$m]')
plt.ylabel(r'$j$ [1/s/cm$^2$]')
plt.title(r'Flux/current')
plt.plot(C7x_microns, C7j)
plt.show()
## Set labels.
plt.xlabel(r'z [$\mu$m]')
plt.ylabel(r'$q$ [W/cm$^2$]')
plt.title(r'Heat flux')
plt.plot(C7x_microns, C7q)
plt.show()

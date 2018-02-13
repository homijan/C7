import numpy as np
from math import pi
from math import exp

## 
kB = 1.3807e-16
me = 9.1094e-28

ne = 5.0e19
Te = 10000.0
gradTe = -1.0
Zbar = 4.0
sigma = 1e15

import argparse
## Create parser object.
parser = argparse.ArgumentParser(description='Analyze diffusive asymptotic of AWBS + compare to C7.')
## Define input arguments.
parser.add_argument("-n", "--ne", help="Electron density at the point.", type=float)
parser.add_argument("-T", "--Te", help="Temperature at the point.", type=float)
parser.add_argument("-g", "--gradTe", help="Temperature gradient at the point.", type=float)
parser.add_argument("-s", "--sigma", help="Electro-ion cross-section.", type=float)
parser.add_argument("-Z", "--Zbar", help="Ionization at the point.", type=float)

## Parse arguments.
args = parser.parse_args()
if args.ne:
    ne = args.ne
if args.Te:
    Te = args.Te
if args.gradTe:
    gradTe = args.gradTe
if args.sigma:
    sigma = args.sigma
if args.Zbar:
    Zbar = args.Zbar

## We use ion density as reference.
ni = ne / Zbar

def vTh(T): 
    return (kB*T/me)**0.5

def fM(v, T):
    return ne/(vTh(T)**3.0*(2.0*pi)**1.5)*exp(-v**2.0/2.0/vTh(T)**2.0)

def rhs(v, T, gradT, Z, E):
    return Z/v*((v**2.0/2.0/vTh(T)**2.0 - 1.5)*gradT/T - E/vTh(T)**2.0)*fM(v, T)

def func(y, v, T, gradT, Z, E):
    f1, = y
    dydt = [(Z - 4.0)/v*f1 + rhs(v, T, gradT, Z, E)]
    return dydt

def solve_fweuler(v, f0, T, gradT, Z, E):
    N = len(v)
    f1 = np.zeros(N)
    f1[0] = f0
    for i in range(N-1):
        dv = v[i+1] - v[i]
        vp = v[i]
        #f1[i+1] = f1[i] + dv*rhs(vp, T, gradT, Z, E)
        f1[i+1] = (1.0 + dv*(Z - 4.0)/vp)*f1[i] + dv*rhs(vp, T, gradT, Z, E)
    return f1

def solve_bweuler(v, f0, T, gradT, Z, E):
    N = len(v)
    f1 = np.zeros(N)
    f1[0] = f0
    for i in range(N-1):
        dv = v[i+1] - v[i]
        vp = v[i]
        f1[i+1] = (f1[i] + dv*rhs(vp, T, gradT, Z, E))/(1.0 + dv*(4.0 - Z)/vp)
    return f1

ml_max = 10.0
ml_min = 0.05
## The heat flux after integration takes the form
## qH = me/Zbar/sigma*128/(2*pi)**0.5*(kB/me)**(7/2)*T**(5/2)*gradT,
## where mfp_ei = v**4/sigma/ni/Zbar, i.e. sigma corresponds to ei collisions.
corr = (688.9*Zbar + 114.4)/(Zbar**2.0 + 1038*Zbar + 474.1)
print "Zbar, corr:", Zbar, corr
cmag = 1./corr
## Classical Lorentz approximation electric field.
Efield = vTh(Te)**2.0*2.5*gradTe/Te

## Solve the AWBS ODE problem.
Nexpl = 5000
vexpl = np.linspace(ml_max*vTh(Te), ml_min*vTh(Te), Nexpl)
dvexpl = (ml_max - ml_min)*vTh(Te)/(Nexpl-1)
Nimpl = 5000
vimpl = np.linspace(ml_max*vTh(Te), ml_min*vTh(Te), Nimpl)
dvimpl = (ml_max - ml_min)*vTh(Te)/(Nimpl-1)
## Use both explicit and implicit method of ODE solve.
sol_expl = solve_fweuler(vexpl, 0.0, Te, gradTe, Zbar, Efield)
sol_impl = solve_bweuler(vimpl, 0.0, Te, gradTe, cmag*Zbar, Efield)

## Post-process transport values
SHf1 = np.zeros(Nimpl)
SHj = np.zeros(Nimpl)
SHq = np.zeros(Nimpl)
SHJ = 0.0
SHQ = 0.0
AWBSf1_impl = np.zeros(Nimpl)
AWBSj_impl = np.zeros(Nimpl)
AWBSq_impl = np.zeros(Nimpl)
AWBSJ_impl = 0.0
AWBSQ_impl = 0.0
for i in range(Nimpl):
    vp = vimpl[i]
    ## The mean free path has standard v^4 dependence, sigma is cross section
    ## given by model and Zbar increases the effect of Coulomb potential in
    ## ei collisions.
    mfp_ei = vp**4.0/sigma/ni/Zbar
    dv = dvimpl
    SHf1[i] = - (Zbar + 0.24)/(Zbar + 4.2)*((vp**2.0/2.0/vTh(Te)**2.0 - 1.5)*gradTe/Te - Efield/vTh(Te)**2.0)*fM(vp, Te)*vp*vp
    SHj[i] = mfp_ei*vp*SHf1[i]
    SHq[i] = mfp_ei*me/2.0*vp*vp*vp*SHf1[i]
    SHJ = SHJ + 4.0*pi/3.0*SHj[i]*dv
    SHQ = SHQ + 4.0*pi/3.0*SHq[i]*dv
    AWBSf1_impl[i] = sol_impl[i]*vp*vp
    AWBSj_impl[i] = mfp_ei*vp*AWBSf1_impl[i]
    AWBSq_impl[i] = mfp_ei*vp*vp*vp*me/2.0*AWBSf1_impl[i]
    AWBSJ_impl = AWBSJ_impl + 4.0*pi/3.0*AWBSj_impl[i]*dv
    AWBSQ_impl = AWBSQ_impl + 4.0*pi/3.0*AWBSq_impl[i]*dv
AWBSf1_expl = np.zeros(Nexpl)
AWBSj_expl = np.zeros(Nexpl)
AWBSq_expl = np.zeros(Nexpl)
AWBSJ_expl = 0.0
AWBSQ_expl = 0.0
for i in range(Nexpl):
    vp = vexpl[i]
    mfp_ei = vp**4.0/sigma/ni/Zbar
    dv = dvexpl
    AWBSf1_expl[i] = sol_expl[i]*vp*vp
    AWBSj_expl[i] = mfp_ei*vp*AWBSf1_expl[i]
    AWBSq_expl[i] = mfp_ei*vp*vp*vp*me/2.0*AWBSf1_expl[i]
    AWBSJ_expl = AWBSJ_expl + 4.0*pi/3.0*AWBSj_expl[i]*dv
    AWBSQ_expl = AWBSQ_expl + 4.0*pi/3.0*AWBSq_expl[i]*dv

## Compare with C7 results
## No explicit treatment of Efield, we use mimicing by reducing ne in source.
C7v, C7mehalff1v5 = np.loadtxt('fe_point_Emimic.txt',  usecols=(0, 4), unpack=True)
C7Q = 0.0
NC7 = C7v.size - 1
for i in range(NC7):
    dC7v = C7v[i] - C7v[i+1]
    C7Q = C7Q + C7mehalff1v5[i]*dC7v
## Explicit treatment of Efield.
C7Ev, C7Emehalff1v5 = np.loadtxt('fe_point_Ecorrect.txt',  usecols=(0, 4), unpack=True)
C7EQ = 0.0
NC7E = C7Ev.size - 1
for i in range(NC7E):
    dC7Ev = C7Ev[i] - C7Ev[i+1]
    C7EQ = C7EQ + C7Emehalff1v5[i]*dC7Ev

## Analytical formula from AWBShos.pdf, providing the Lorentz gas approximation
## further multiplied by SH low Z factor.
mfp_ei = (vTh(Te))**4.0/sigma/ni/Zbar
L = Te / abs(gradTe)
SHQ_analytic = - (Zbar + 0.24)/(Zbar + 4.2) * 128.0/(2.0*pi)**0.5*ne*vTh(Te)*kB*Te*mfp_ei*gradTe/Te
Kn = mfp_ei/L
Kn_flux = SHQ_analytic / ((Zbar + 0.24)/(Zbar + 4.2) * 128.0/(2.0*pi)**0.5 * ne * vTh(Te) * kB * Te)
## Show the Knudsen number
print 'Kn: ', Kn, 'Kn from flux: ', Kn_flux 
## Express flux proportionality with respect to SHQ_analytic
proporC7EQ = C7EQ / SHQ_analytic
proporC7Q = C7Q / SHQ_analytic

## Print integrated values
print "SHQ:          ", SHQ
print "SHQ_analytic: ", SHQ_analytic
print "C7EQ:         ", C7EQ
print "C7Q:          ", C7Q
print "AWBSQ_impl:   ", AWBSQ_impl
print "AWBSQ_expl:   ", AWBSQ_expl
print "SHJ:        ", SHJ
print "AWBSJ_impl: ", AWBSJ_impl
print "AWBSJ_expl: ", AWBSJ_expl

# Physical fix by a magic constant.
# Z = 1
#cmag = 1.882
# Z = 5
#cmag = 1.59575
# Z = 10
#cmag = 1.565
# Z = 20
#cmag = 1.56
# Z = 50
#cmag = 1.585
# Z = 100
#cmag = 1.65
# Z = 200
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
#print icoeffs
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

import matplotlib.pyplot as plt
import matplotlib
font = {'family' : 'Sans',
        #'weight' : 'bold',
        'size'   : 16}
matplotlib.rc('font', **font)

#plt.plot(Zbars, coeffs, 'rx', label='pointwise corrections')
#plt.plot(ZZ, icoeffs, 'b', label=str(p1)+r'$Z + $'+str(p2)+r'$/Z^2 + $'+str(q1)+r'$Z + $'+str(q2))
#plt.legend(loc='best')
#plt.xlabel('Z')
#plt.title('Rational function fit')

## Shrink the x axis appropriately.
mult = 8
p_vimpl = vimpl[vimpl < mult*vTh(Te)]
p_SHq = SHq[vimpl < mult*vTh(Te)]
p_AWBSq_impl = AWBSq_impl[vimpl < mult*vTh(Te)]
p_vexpl = vexpl[vexpl < mult*vTh(Te)]
p_AWBSq_expl = AWBSq_expl[vexpl < mult*vTh(Te)]

## Set labels.
plt.ylabel(r'$q_1 = m_e/2 v^2 f_1 v^2$ [a.u.]')
plt.xlabel('v/vT')
plt.title('Z = '+str(Zbar)+', Kn = '+"{:.1e}".format(Kn))
## Plot kinetic analysis.
plt.plot(p_vimpl/vTh(Te), p_SHq, 'b', label='SH')
#plt.plot(p_vexpl/vTh(Te), p_AWBSq_expl, 'g-.', label='AWBS')
plt.plot(p_vimpl/vTh(Te), p_AWBSq_impl, 'r--', label=r'AWBS$^{*}$')
if (len(C7Ev)<30):
   plt.plot(C7Ev/vTh(Te), C7Emehalff1v5 / (4.0*pi/3.0), 'kx', label='C7E('+"{:.2f}".format(proporC7EQ)+r'q$_{SH}$)')
else:
   plt.plot(C7Ev/vTh(Te), C7Emehalff1v5 / (4.0*pi/3.0), 'k', label='C7E('+"{:.2f}".format(proporC7EQ)+r'q$_{SH}$)')
#plt.plot(C7v/vTh(Te), 1.5 * C7mehalff1v5 / (4.0*pi/3.0), 'k:', label=r'C7')
plt.plot(C7v/vTh(Te), C7mehalff1v5 / (4.0*pi/3.0), 'k--', label=r'C7$^{*}$('+"{:.2f}".format(proporC7Q)+r'q$_{SH}$)')
plt.legend(loc='best')
#plt.grid()
plt.show()

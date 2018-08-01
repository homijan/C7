## Input plasma parameters in cgs (eV).
ne = 1e19
Te = 25.0
Zbar = 3.5
beta = 4.0
B_tesla = 15.0
A = 3.5

## Collisional plasma.
#coulLog = 2.5
## Laser corona plasma.
coulLog = 7.09
## Low density astro-plasma.
#coulLog = 10.0

###############################################################################
########### Essential collisional physics quantities in cgs ###################
###############################################################################
## Fundamental physical constants in cgs units. 
kB = 1.6022e-12
me = 9.1094e-28
mp = 1.6726e-24
qe = 4.8032e-10
c = 2.9979e10
## Classical v^4 cross-section in cgs units.
crs = 8.1027575e17 ## Matching the SH diffusive flux.
## Electron thermal velocity.
def vTe(T): 
    return (kB*T/me)**0.5
## Ion thermal velocity.
def vTi(A, T): 
    return (kB*T/A/mp)**0.5
## Maxwell-Boltzmann distribution.
def fM(ne, v, T):
    return ne/(vTe(T)**3.0*(2.0*pi)**1.5)*exp(-v**2.0/2.0/vTe(T)**2.0)
###############################################################################
###############################################################################

## Test particle velocity anc classical collision model.
v_electron = beta * vTe(Te) 
#nu_tot = crs * coulLog * ne * (Zbar+1.0) / v_electron**3.0
nu_LMV = crs * coulLog * ne * (Zbar+1.0)**0.5 / v_electron**3.0
#mfp_tot = (vTe(Te))**4.0/crs/coulLog/ne/(Zbar+1.0)
#mfp_LMV = v_electron**4.0 / crs / coulLog / ne / (Zbar+1.0)**0.5
mfp = v_electron / nu_LMV

## The same classical collision model expressed in plasma parameters.
#C_nu = crs * me**1.5 / kB**1.5
C_nu = 1.09848064563e-05
nu_ = C_nu * ne * (Zbar+1.0)**0.5 * coulLog / (beta**3.0 * Te**1.5) 
#C_mfp = kB**2.0 / me**2.0 * 1. / crs
C_mfp = 3.81786939841e+12
mfp_ = C_mfp * beta**4.0 * Te**2.0 / (ne * (Zbar + 1.0)**0.5 * coulLog)

#f_B = 2.8e10 * B_tesla # wikipedia :)
#omega_B = qe / me / c * B_tesla * 1e4 
#f_B_ = omega_B / 2.0 / 3.14159
#radius = v_electron / omega_B
#radius = (kB / me)**0.5 * me * c / qe * beta * Te**0.5 / B_tesla / 1e4
#C_r = (kB / me)**0.5 * me * c / qe / 1e4
C_r = 2.38445824407e-4
radius_ = C_r * beta * Te**0.5 / B_tesla

## Test particle velocity anc classical collision model.
v_ion = vTi(A, Te) 
nu_ii = crs * me**2.0 / (A*mp)**2.0 * coulLog * ne * Zbar**3.0 / v_ion**3.0
mfp_ii = v_ion / nu_ii

#print "nu: ", nu_LMV, nu_
#print "mfp[microns]: ", mfp * 1e4, mfp_ * 1e4 
#print "C_nu, C_mfp: ", C_nu, C_mfp
#print "f_B: ", f_B, f_B_ 

print "ne[cm-3]: ", ne
print "Te[eV]:   ", Te
print "Zbar:     ", Zbar
print "beta:     ", beta, "...electron thermal velocity multiple."
print "B[tesla]: ", B_tesla
print "A:        ", A
print "mfp[microns] = C_mfp * beta^4 * Te^2 / ne / (Z+1)^0.5: "
print mfp_ * 1e4
print "radius[microns] = C_r * beta * Te^0.5 / B: "
print radius_ * 1e4
print "C_mfp, C_r: ", C_mfp * 1e4, C_r * 1e4
print ""
print "Electron velocity[km/s] (beta * thermal velocity):", beta * vTe(Te) / 1e4
print "Ion thermal velocity[km/s]:", vTi(A, Te) / 1e4
print ""
print "Thermal mfpe_ii[microns]: ", mfp_ii * 1e4
print "shock thickness[microns]: ", mfp_ii * (A * mp / me)**0.5 * 1e4

import numpy as np
from math import pi
from matplotlib import pyplot as plt

## Pick a spatial position index equals microns.
#point = 500 # reference for scaling
#point = 685
point = 797
Zbar = 2.0

## Fundamental physical constants in cgs. 
kB = 1.6022e-12
me = 9.1094e-28
qe = 4.8032e-10
c = 2.99792458e10
## Cross-section attributes.
sigma = 8.1027575e17
coulLog = 7.09

## You need to know the data-axis size and range.
Nv = 500
Nx = 1000
x_min = 0.0 * 1e-4 ## cm
x_max = 5902.4 / 2.0 / pi * 1e-4 ## cm
gamma_min = 1.0
gamma_max = 1.04
## Scale velocity and spatial axis.
x = np.linspace(x_min, x_max, Nx)
gamma = np.linspace(gamma_min, gamma_max, Nv)
v = np.array((1.0 - 1.0 / gamma**2.0)**0.5 * c)
## Load and transform the data.
## Calder corresponding temperature in eV.
_Te = np.loadtxt('Te/Te_00110000.txt') * 1e3
## Distribution function f0 and f1.
f0 = np.zeros((Nv, Nx))
f1 = np.zeros((Nv, Nx))
q_data = np.loadtxt('f/qxnr_00110000.txt')
a_data = np.loadtxt('f/axan_00110000.txt')
## Store the distribution function on a grid.
for i in range(Nv):
    for j in range(Nx):
        f0[i, j] = q_data[Nx * i + j]
        f1[i, j] = a_data[Nx * i + j]
## A given point distribution.
_f0_point = f0[:, point]
_f1_point = f1[:, point]

## PIC data are noisy, so some temperature smoothing is necessary.
from scipy.interpolate import splev, splrep 
## Obtain smooth profiles of temperature and its gradient.
smooth = 700 # lower equals less smoothing
## Find a spline for the temperature data.
Tetck = splrep(x, _Te, s=smooth)
## Assign smooth profiles profiles.
Te = splev(x, Tetck, der=0)
gradTe = splev(x, Tetck, der=1)
#print "x_max: ", x_max
#print "max(_Te), max(Te): ", max(_Te), max(Te)
## Obtain profiles of distribution.
smooth = 0 # lower equals less smoothing
## Find a spline for the distribution data.
f0tck = splrep(v, _f0_point, s=smooth)
f1tck = splrep(v, _f1_point, s=smooth)
## Assign smooth profiles profiles.
f0_point = splev(v, f0tck, der=0) * v**2.0
f1_point = splev(v, f1tck, der=0) * v**2.0
j1_point = me / 2.0 * f1_point * v
q1_point = me / 2.0 * f1_point * v**3.0

## Thermal velocity.
vTe = (kB * Te / me)**0.5
vTe_point = vTe[point]
## SH flux.
## Zbar correction.
SHcorr = (Zbar + 0.24)/(Zbar + 4.2)
SHQ_Wcm2_cgs_eV = - SHcorr * 1.31e3 / coulLog / Zbar * Te**2.5 * gradTe

## Plot results.
#plt.plot(x, _Te)
#plt.plot(x, Te)
#plt.show()
#plt.plot(x, gradTe)
#plt.show()
plt.plot(x, SHQ_Wcm2_cgs_eV)
plt.show()
x_point = x[point] * 1e4
Te_point = Te[point]
print "x(point)[microns]: ", x_point
print "Te(point)[eV]: ", Te[point]
dv = v[1] - v[0]
ne_point = sum(f0_point) * dv
## Valid explicitly for simulation run ne = 5e20 cm-3.
scale = 4.62158960888e-11
#scale = 5. / 7. * 1e-10
#scale = 5e20 / ne_point
print "scale: ", scale
print "ne(point): ", ne_point * scale
## Scale EDF.
f0_point = f0_point * scale
f1_point = f1_point * scale
j1_point = j1_point * scale
q1_point = q1_point * scale


plt.plot(v / vTe_point, f0_point)
plt.show()
plt.plot(v / vTe_point, f1_point)
plt.show()
plt.plot(v / vTe_point, j1_point)
plt.show()
q_point = sum(q1_point) * dv
print "q(point)[Wcm2]: ", q_point * 1e-7

import scipy.fftpack as fftp
def fourier_series(x, y, wn, n=None):
    # get FFT
    myfft = fftp.fft(y, n)
    # kill higher freqs above wavenumber wn
    myfft[wn:-wn] = 0
    # make new series
    y2 = fftp.ifft(myfft)
    plt.figure(num=None)
    plt.plot(x, y, x, y2.real)
    plt.show()
    return y2.real

y = q1_point
wn = 30
q1_point = fourier_series(v, y, wn)
#q1_point[v / vTe_point > 11] = 0.0
plt.plot(v / vTe_point, q1_point)
plt.show()

np.savetxt('f/F0F1x_Calder_Zeq2_tanh_50mic_f0f1_1.10e-11_'+'{:.1f}'.format(x_point)+'mic.txt', np.transpose([v, f0_point, q1_point]))

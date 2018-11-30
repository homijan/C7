import numpy as np
## Graphics.
import matplotlib.pyplot as plt
import matplotlib
## Interpolation.
from scipy.interpolate import splev, splrep

## LTE functions as fM.
from LTE import *

def getdata(filename):
    data = np.array(np.loadtxt(filename))
    x = data[:, 0]
    y = data[:, 1]
    return x, y

## Auxiliary function operating on derivatives.
def map_f(xs, fs, xmap):
    ## Map the (xs, fs) interpolation on the xmap points.
    smooth = 0 # lower less smoothing
    ## Find a spline for the f data.
    ftck = splrep(xs, fs, s=smooth)
    return splev(xmap, ftck, der=0)

### Z=1
filename_Cal = 'Calder_Tenoisy_Z1_ne5e20_tanh50_20ps.txt'
x_Cal, T_Cal = getdata(filename_Cal)
## A silly change avoiding SNB bad resultsi :).
x_Cal[0] = 0.0
## PIC data are noisy, so some temperature smoothing is necessary.
T_noisy = T_Cal
## Obtain smooth profiles of temperature and its gradient.
smooth = 0.0001 # lower equals less smoothing
## Find a spline for the temperature data.
Tetck = splrep(x_Cal, T_Cal, s=smooth)
## Assign smooth profiles profiles.
T_Cal = splev(x_Cal, Tetck, der=0)
## Save the smoothed (only high frequency low amplitude noise).
np.savetxt('Calder_Te_Z1_ne5e20_tanh50_20ps.txt', np.transpose([x_Cal, T_Cal]))

Te_460 = 1e3 * map_f(x_Cal, T_Cal, 460.0)
Te_580 = 1e3 * map_f(x_Cal, T_Cal, 580.0)
vT_460 = vTh(Te_460)
vT_580 = vTh(Te_580)
print("Te_460 [eV]:", Te_460)
print("Te_580 [eV]:", Te_580)

plt.plot(x_Cal, T_Cal, label='T-Calder')
plt.plot(x_Cal, T_noisy, label='T-Calder-noisy')
plt.legend()
plt.show()

filename_Cal = 'Calder_Qx_Z1_ne5e20_tanh50_20ps.txt'
x_Cal, Q_Cal = getdata(filename_Cal)
Q_460 = map_f(x_Cal, Q_Cal, 460.0)
Q_580 = map_f(x_Cal, Q_Cal, 580.0)

plt.plot(x_Cal, Q_Cal, label='Q-Calder')
plt.legend()
plt.show()

filename_Cal = 'Calder_f1z1_20ps_460mu.txt'
x_Cal, f1_Cal = getdata(filename_Cal)
## Proper cgs units.
intOmega = 3.0 / 4.0 / np.pi #4.0 * np.pi / 3.
v_Cal = x_Cal * vT_460
dv = v_Cal[1] - v_Cal[0]
q_ = dv * sum(f1_Cal)
#scale = 1e7 * Q_460 / intOmega / q_
# unique value for this current data.
scale = 677120186180373.6
f1_Cal *= scale #vT_460**5.0
## Check the f1.
Q_ = intOmega * dv * sum(f1_Cal)
print("Q_460 [erg/s/cm2], Q_:", 1e7 * Q_460, Q_)
print("scale:", scale)

np.savetxt('Calder_f1_Z1_ne5e20_tanh50_20ps_460mu.txt', np.transpose([v_Cal, f1_Cal]))

plt.plot(x_Cal, f1_Cal, label='f1-460-Calder')
plt.legend()
plt.show()

filename_Cal = 'Calder_f1z1_20ps_580mu.txt'
x_Cal, f1_Cal = getdata(filename_Cal)
## Proper cgs units.
v_Cal = x_Cal * vT_580
dv = v_Cal[1] - v_Cal[0]
q_ = dv * sum(f1_Cal)
#scale = 1e7 * Q_580 / intOmega / q_
# unique value for this current data.
scale = 677120186180373.6
f1_Cal *= scale #me * vT_580**5.0
## Check the f1.
Q_ = intOmega * dv * sum(f1_Cal)
print("Q_580 [erg/s/cm2], Q_:", 1e7 * Q_580, Q_)
print("scale:", scale)

np.savetxt('Calder_f1_Z1_ne5e20_tanh50_20ps_580mu.txt', np.transpose([v_Cal, f1_Cal]))

plt.plot(x_Cal, f1_Cal, label='f1-580-Calder')
plt.legend()
plt.show()

### Z=10
filename_Cal = 'Calder_Te_Z10_ne5e20_tanh50_12ps.txt'
x_Cal, T_Cal = getdata(filename_Cal)
Te_460 = 1e3 * map_f(x_Cal, T_Cal, 460.0)
Te_580 = 1e3 * map_f(x_Cal, T_Cal, 580.0)
vT_460 = vTh(Te_460)
vT_580 = vTh(Te_580)
print("Te_460 [eV]:", Te_460)
print("Te_580 [eV]:", Te_580)

plt.plot(x_Cal, T_Cal, label='T-Calder')
plt.legend()
plt.show()

filename_Cal = 'Calder_Qx_Z10_ne5e20_tanh50_12ps.txt'
x_Cal, Q_Cal = getdata(filename_Cal)
Q_460 = map_f(x_Cal, Q_Cal, 460.0)
Q_580 = map_f(x_Cal, Q_Cal, 580.0)

plt.plot(x_Cal, Q_Cal, label='Q-Calder')
plt.legend()
plt.show()

filename_Cal = 'Calder_f1z10_12ps_460mu.txt'
x_Cal, f1_Cal = getdata(filename_Cal)
## Proper cgs units.
intOmega = 3.0 / 4.0 / np.pi #4.0 * np.pi / 3.
v_Cal = x_Cal * vT_460
dv = v_Cal[1] - v_Cal[0]
q_ = dv * sum(f1_Cal)
scale = 1e7 * Q_460 / intOmega / q_
f1_Cal *= scale #me * vT_460**5.0
## Check the f1.
Q_ = intOmega * dv * sum(f1_Cal)
print("Q_460 [erg/s/cm2], Q_:", 1e7 * Q_460, Q_)
print("scale:", scale)

np.savetxt('Calder_f1_Z10_ne5e20_tanh50_12ps_460mu.txt', np.transpose([v_Cal, f1_Cal]))

plt.plot(x_Cal, f1_Cal, label='f1-460-Calder')
plt.legend()
plt.show()

filename_Cal = 'Calder_f1z10_12ps_580mu.txt'
x_Cal, f1_Cal = getdata(filename_Cal)
## Proper cgs units.
v_Cal = x_Cal * vT_580
dv = v_Cal[1] - v_Cal[0]
q_ = dv * sum(f1_Cal)
scale = 1e7 * Q_580 / intOmega / q_
f1_Cal *= scale #me * vT_580**5.0
## Check the f1.
Q_ = intOmega * dv * sum(f1_Cal)
print("Q_580 [erg/s/cm2], Q_:", 1e7 * Q_580, Q_)
print("scale:", scale)

np.savetxt('Calder_f1_Z10_ne5e20_tanh50_12ps_580mu.txt', np.transpose([v_Cal, f1_Cal]))

plt.plot(x_Cal, f1_Cal, label='f1-580-Calder')
plt.legend()
plt.show()

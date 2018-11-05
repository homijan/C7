import numpy as np

## Fundamental physical constants in cgs. 
kB = 1.6022e-12
me = 9.1094e-28
qe = -4.8032e-10

## Basic local thermal equilibrium functions.
def vTh(T): 
    return (kB*T/me)**0.5
def fM(v, n, T):
    return n/(vTh(T)**3.0*(2.0*np.pi)**1.5)*np.exp(-v**2.0/2.0/vTh(T)**2.0)
def dfMdv(v, n, T):
    return - v/vTh(T)**2.0*fM(v, n, T)
def dfMdx(v, n, T, dndx, dTdx):
    dfMdn = fM(v, n, T) / n
    dfMdT = (0.5*v**2.0 / vTh(T)**2.0 - 3.0 / 2.0) / T * fM(n, T, v)
    return dfMdn * dndx + dfMdT * dTdx

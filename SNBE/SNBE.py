import numpy as np
## Interpolation.
from scipy.interpolate import splev, splrep
## Brent solver for corr factors.
from scipy import optimize
## delta corr fit.
from scipy.optimize import curve_fit
## Graphics.
import matplotlib.pyplot as plt
import matplotlib

## LTE functions as fM.
from LTE import *
## Diffusion solver.
from diffusion_solver import *

## Cross section in cgs.
sigma = 8.1027575e17 ## Matching the SH diffusive flux.
## Coulomb logarithm.
coulLog = 7.09
## Gamma factor cross section.
Gamma = sigma * coulLog
## BGK e-e collision scaling coefficient.
rr = 0.5
# Scaling from microns to cm.
mictocm = 1e-4

## Auxiliary function operating on derivatives.
def map_f(xs, fs, xmap):
    ## Map the (xs, fs) interpolation on the xmap points.
    smooth = 0 # lower less smoothing
    ## Find a spline for the f data.
    ftck = splrep(xs, fs, s=smooth)
    return splev(xmap, ftck, der=0)
def dfdx(xs, fs):
    ## Compute the first derivative from a data interpolation.
    smooth = 0 # lower less smoothing
    ## Find a spline for the f data.
    ftck = splrep(xs, fs, s=smooth)
    dfdx = splev(xs, ftck, der=1)
    #d2fdx2 = splev(xs, ftck, der=2)
    return dfdx

## Collisions and local approximations.
# Auxiliary functions.
def f1_grad(xs, v, n, T, E):
    grad_n = dfdx(xs, n)
    grad_T = dfdx(xs, T)
    return grad_n / n + (v**2.0 / 2.0 / vTh(T)**2.0 - 3.0 / 2.0) * grad_T / T - qe / me / vTh(T)**2.0 * E
def f1_grad_simple(xs, v, n, T):
    grad_n = dfdx(xs, n)
    grad_T = dfdx(xs, T)
    #return grad_n / n + grad_T / T
    return grad_T / T
# Electron-electron collisions.
def nue(v, n):
    return Gamma * n / v**3.0
# Electron-ion collisions.
def nuei(v, n, Z):
    return Z * nue(v, n)
# Spitzer-Harm ionization scaling of heat flux.
def xi(Z):
    return (Z + 0.24) / (Z + 4.2)
# Anisotropic part of EDF in local approximation.
def func_f1M(xs, v, n, T, Z, E):
    return - xi(Z) * v / nuei(v, n, Z) * fM(v, n, T) * f1_grad(xs, v, n, T, E)
def func_df1Mdv(xs, v, n, T, Z, E):
    return - xi(Z) * (lambdaei(v, n, Z, E) * fM(v, n, T) * v / vTh(T)**2.0 * dfdx(xs, T) / T + (dlambdaeidv(v, n, Z, E) * fM(v, n, T) + lambdaei(v, n, Z, E) * dfMdv(v, n, T)) * f1_grad(xs, v, n, T, E))
# Spitzer-Harm heat flux used in hydro codes.
def SH_flux(xs, T, Z):
    return - xi(Z) * 1.31e10 / coulLog / Z * T**2.5 * dfdx(xs, T)

## SNB e-e and e-i mean free paths.
def lambdae(v, n, Z):
    #return 2.0 * v / nue(v, n)
    return xi(Z) * v / rr / nue(v, n) 
def lambdaei(v, n, Z, E):
    mfpei = xi(Z) * v / nuei(v, n, Z)
    #SIMPLE
    if (SIMPLE):
        mfpei = 1.0 / (nuei(v, n, Z) / xi(Z) / v + abs(qe * E) / 0.5 / me / v**2)
    return mfpei
def dlambdaeidv(v, n, Z, E):
    return xi(Z) * 4.0 * v**3.0 / Gamma / n 

## SNBE scheme coefficients.
# Diffusion solver returns solution of d[Ddfdx]dx - A * f = S
def coeff_D(v, n, Z, E):
    return lambdaei(v, n, Z, E) / 3.0
def coeff_A(v, n, Z):
    return 1.0 / lambdae(v, n, Z)
def coeff_S(xs, v, n, T, Z, E):
    grad = f1_grad(xs, v, n, T, E)
    E_f1M = qe / me / v / 3.0 * E * (func_df1Mdv(xs, v, n, T, Z, E) + 2.0 / v * func_f1M(xs, v, n, T, Z, E))
    #SIMPLE
    if (SIMPLE):
        grad = f1_grad_simple(xs, v, n, T)
        E_f1M[:] = 0.0
    f1 = xi(Z) * v / nuei(v, n, Z) / 3.0 * fM(v, n, T) * grad
    div_f1 = dfdx(xs, f1) 
    # Revert sign due to the diffusion solver.
    return - div_f1 + E_f1M

# If run as main program
if __name__ == "__main__":
    # Switch for a simple SNB (E_L) run.
    SIMPLE = True
    #SIMPLE = False
    # Number of x-cells.
    N = 400
    # Number of energy groups "v-cells".
    Ngr = 250
    # Reference values for density and ionization.    
    ne_ref = 5e20
    Zbar_ref = 10.0
    # A reference point for the distribution function output.
    x_point = 580.0 * mictocm
    # Number of consistent electric field iterations.
    E_iter_max = 30

    import argparse
    ps = argparse.ArgumentParser( description = 'SNBE - steady state kinetic solution on given Te, ne, Zbar plasma profiles in [cgs] and T [eV].')
    ps.add_argument( '-N', '--Ncells', type = int, help = 'Number of cells in spatial discretization.' )
    ps.add_argument( '-Ngr', '--Ngroups', type = int, help = 'Number of groups in velocity discretization.' )
    ps.add_argument( '-Tinf', '--Te_inputfile', type = str, help = 'Temperature input file.' )
    ps.add_argument( '-ninf', '--ne_inputfile', type = str, help = 'Electron density input file.' )
    ps.add_argument( '-Zinf', '--Zbar_inputfile', type = str, help = 'Mean ionization input file.' )
    ps.add_argument( '-Qo', '--Q_outname', type = str, help = 'Heat flux Q output file name.' )
    ps.add_argument( '-F1o', '--F1_outname', type = str, help = 'F1 output file name.' )
    ps.add_argument( '-ne', '--ne_value', type = float, help = 'Electron density value.' )
    ps.add_argument( '-Z', '--Zbar_value', type = float, help = 'Mean ionization value.' )
    ps.add_argument( '-pt', '--EDF_outpoint', type = float, help = 'Output point of EDF.' )
    ## A no value argument solution.
    ps.add_argument("-s", "--plotshow", action='store_true', help="Provide plots of kinetic simulation related quantities.")
    # Read args.
    args = ps.parse_args()
    # Reflect appropriate arguments.
    if (args.Ncells):
        N = args.Ncells
    if (args.Ngroups):
        Ngr = args.Ngroups
    if (args.ne_value):
        ne_ref = args.ne_value
    if (args.Zbar_value):
        Zbar_ref = args.Zbar_value
    if (args.EDF_outpoint):
        x_point = args.EDF_outpoint

    # Spatial range.
    x_min = 0.0
    x_max = 700.0 * mictocm
    
    ## If input files provided.
    # Computational domain is set according to the Te input file.
    if (args.Te_inputfile):
        #filename = 'temperature.dat'
        data = np.array(np.loadtxt(args.Te_inputfile))
        xx_TT = data[:, 0]
        TT = data[:, 1]
        x_min = min(xx_TT)
        x_max = max(xx_TT) 
    if (args.ne_inputfile):
        #filename = 'ne.dat'
        data = np.array(np.loadtxt(args.ne_inputfile))
        xx_nn = data[:, 0]
        nn = data[:, 1]
    if (args.Zbar_inputfile):
        #filename = 'zbar.dat'
        data = np.array(np.loadtxt(args.Zbar_inputfile))
        xx_ZZbar = data[:, 0]
        ZZbar = data[:, 1]

    ## Staggered scheme.
    # Diffusive coefficient.
    xD = np.linspace(x_min, x_max, N-1)
    # Uniform mesh cell size.
    dx = xD[1] - xD[0]
    # Cell centered values.
    x = np.linspace(x_min - dx/2.0, x_max + dx/2.0, N)
    
    ## Diffusion scheme coefficients.
    D = np.zeros(N-1)
    A = np.zeros(N)
    S = np.zeros(N)
    f1M = np.zeros(N-1)
    ## Distribution function. 
    EDF_df0dx = np.zeros((Ngr, N-1))
    EDF_f1 = np.zeros((Ngr, N-1))
    EDF_f1M = np.zeros((Ngr, N-1))
    EDF_q1M = np.zeros((Ngr, N-1))
    EDF_q1 = np.zeros((Ngr, N-1))
    EDF_j1 = np.zeros((Ngr, N-1))

    ## Plasma profiles.
    # Constant profiles of ionization and electron density.
    Zbar = Zbar_ref * np.ones(N)
    Zbar_xD = Zbar_ref * np.ones(N-1)
    ne = ne_ref * np.ones(N)
    ne_xD = ne_ref * np.ones(N-1)
    # Heat-bath profile of temperature.
    s = 1.0 / 25.0
    Te = 1e3 * (0.575 - 0.425 * np.tanh((x - 450.0 * mictocm) * s / mictocm))
    Te_xD = 1e3 * (0.575 - 0.425 * np.tanh((xD - 450.0 * mictocm) * s / mictocm))
    #plt.plot(xD, Te_xD)
    #plt.show()
    
    # Mapped plasma profiles from input files.
    if (args.Te_inputfile):
        Te = map_f(xx_TT, TT, x)
        Te_xD = map_f(xx_TT, TT, xD)
    if (args.ne_inputfile):
        ne = map_f(xx_nn, nn, x)
        ne_xD = map_f(xx_nn, nn, xD)
    if (args.Zbar_inputfile):
        Zbar = map_f(xx_ZZbar, ZZbar, x)
        Zbar_xD = map_f(xx_ZZbar, ZZbar, xD)


    ## Test point structures.
    # Index of the x_point.
    index_point = int((x_point - x_min) / dx)
    vTh_point = vTh(Te_xD[index_point])
 
    ## Velocity space range.
    v_max_global = 7.0 * vTh(max(Te))
    # Equidistant stepping.
    dv = v_max_global / Ngr
    # Cell centered min and max velocity.
    v_min = dv / 2.0
    v_max = v_max_global - dv / 2.0
    # v-space discretization.
    v_grs = np.linspace(v_min, v_max, Ngr) 

    ## Electric field structures.
    # Reference Lorentz electric field.
    E_L = me / qe * vTh(Te)**2.0 * (dfdx(x, ne) / ne + 2.5 * dfdx(x, Te) / Te)
    E_L_xD = me / qe * vTh(Te_xD)**2.0 * (dfdx(xD, ne_xD) / ne_xD + 2.5 * dfdx(xD, Te_xD) / Te_xD)  
    ## Constant denominator in E field definition.
    E_denom_xD = np.zeros(N-1)
    # Nonlocal perturbation of E field.
    deltaE_xD = np.zeros(N-1)

    ## Initial electric field given by Lorentz E field.
    E_xD = E_L_xD

    ## Solve SNB with E field relaxation.
    if (SIMPLE):
        E_iter_max = 1
    print "Running SNBE on", N, "cells and", Ngr,"groups."
    for iter in range(E_iter_max):
        ## Calculate the electric field.   
        # Constant denominator in the E field definition.
        E_denom_xD[:] = 0.0
        for gr in range(Ngr): 
            v = v_grs[gr]
            E_denom_xD[:] += qe / me * lambdaei(v, ne_xD, Zbar_xD, E_xD) * dfMdv(v, ne_xD, Te_xD) * v**2.0 * dv
        # Nonlocal contribution in the E field definition.
        deltaE_xD[:] = 0.0
        for gr in range(Ngr): 
            v = v_grs[gr]
            deltaE_xD[:] += - v * lambdaei(v, ne_xD, Zbar_xD, E_xD) * EDF_df0dx[gr][:] * v**2.0 * dv / E_denom_xD
        # Actual value of the electric field.
        E_xD = E_L_xD + deltaE_xD 
        # Map E_xD on x mesh.
        E = map_f(xD, E_xD, x)

        # Diffusion solver over the entire velocity space. 
        for gr in range(Ngr): 
            v = v_grs[gr]
            D[:] = coeff_D(v, ne_xD, Zbar_xD, E_xD)
            A[:] = coeff_A(v, ne, Zbar)
            S[:] = coeff_S(x, v, ne, Te, Zbar, E)
            f1M[:] = func_f1M(xD, v, ne_xD, Te_xD, Zbar_xD, E_xD) 

            f, df0dx = IsolatedDiffusionProblem(D, A, S, dx)
            # Scaling, since diffusion solver uses D to evaluate q.
            EDF_df0dx[gr][:] = df0dx
            EDF_f1[gr][:] = - D * df0dx * 3.0
            EDF_f1M[gr][:] = f1M[:]
            # Flux moment EDF.
            EDF_q1[gr][:] = me / 2.0 * v**2.0 * v * EDF_f1[gr][:] * v**2.0
            EDF_q1M[gr][:] = me / 2.0 * v**2.0 * v * EDF_f1M[gr][:] * v**2.0
            # Current moment.
            EDF_j1[gr][:] = qe * v * (EDF_f1[gr][:] + EDF_f1M[gr][:]) * v**2.0

    # Heat flux integration of contributions from nonlocal and local EDF parts.
    Qnonlocal = 4.0 * np.pi / 3.0 * np.sum(EDF_q1, 0) * dv
    QM = 4.0 * np.pi / 3.0 * np.sum(EDF_q1M, 0) * dv
    J = 4.0 * np.pi / 3.0 * np.sum(EDF_j1, 0) * dv

    ## Prepare heat flux and EDF for output.
    # Sum of local and nonlocal parts plus unit scaling.
    Q_Wcm2 = (QM + Qnonlocal) * 1e-7
    # EDF at given point.
    q1_point = np.zeros(Ngr)
    q1M_point = np.zeros(Ngr)
    for gr in range(Ngr):
        q1_point[gr] = EDF_q1[gr][index_point]
        q1M_point[gr] = EDF_q1M[gr][index_point]
 
    if (args.plotshow):
        #plt.plot(xD, QM * 1e-7, label='QM')
        plt.plot(xD, SH_flux(xD, Te_xD, Zbar_xD) * 1e-7, label='QSH')
        plt.plot(xD, (QM + Qnonlocal) * 1e-7, label='QM + Q')
        plt.legend()
        plt.show()

        plt.plot(xD, E_L_xD, label=r'$E_L$')
        plt.plot(xD, E_xD, label=r'$E_L + \delta E$')
        plt.legend()
        plt.show()

        plt.plot(xD, J, label='current')
        plt.legend()
        plt.show()

        #plt.plot(v_grs / vTh_point, q1_point, label='q1')
        plt.plot(v_grs / vTh_point, q1M_point, label='q1M')
        plt.plot(v_grs / vTh_point, q1M_point + q1_point, label='q1M + q1')
        plt.legend()
        plt.show()
   
    ## Save output. 
    # Heat flux, x axis in microns.
    if (args.Q_outname):
        np.savetxt(args.Q_outname, np.transpose([xD * 1e4, Q_Wcm2]))
    # EDF.
    if (args.F1_outname):
        ## Find a spline for the data.
        smooth = 0 # lower less smoothing
        q1_tck = splrep(v_grs, q1M_point + q1_point, s=smooth)
        # fine data
        Nfine = 10000
        v_fine = np.linspace(min(v_grs), max(v_grs), Nfine)
        q1_fine = splev(v_fine, q1_tck, der=0)
        ## Store fine data
        np.savetxt(args.F1_outname, np.transpose([v_fine, q1_fine]))
        #plt.plot(v_fine, q1_fine)
        #plt.show()

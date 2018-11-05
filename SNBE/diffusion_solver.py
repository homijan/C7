#!/usr/bin/env python

import numpy as np

## Tri Diagonal Matrix Algorithm(a.k.a Thomas algorithm) solver
def TDMAsolver(a, b, c, d):
    '''
    TDMA solver, a b c d can be NumPy array type or Python list type.
    refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    and to http://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm)
    '''
    nf = len(d) # number of equations
    ac, bc, cc, dc = map(np.array, (a, b, c, d)) # copy arrays
    for it in xrange(1, nf):
        mc = ac[it-1]/bc[it-1]
        bc[it] = bc[it] - mc*cc[it-1] 
        dc[it] = dc[it] - mc*dc[it-1]
        	    
    xc = bc
    xc[-1] = dc[-1]/bc[-1]

    for il in xrange(nf-2, -1, -1):
        xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]

    return xc

def IsolatedDiffusionProblem(D, A, S, dx):
    # Prepare lower-diagonal-upper elements.
    N = len(S)
    a = np.zeros(N-1) 
    b = np.zeros(N)
    c = np.zeros(N-1)
    d = S * dx * dx
    # Diagonal elements.
    for i in range(1, N-1):
        b[i] = - (D[i] + D[i-1] + A[i] * dx * dx)
    # Lower-diagonal elements.
    for i in range(0, N-2):
        a[i] = D[i]
    # Upper-diagonal elements.
    for i in range(1, N-1):
        c[i] = D[i]
    # Neumann BCs.
    b[0] = -1.0
    c[0] = 1.0
    d[0] = 0.0
    a[N-2] = 1.0
    b[N-1] = -1.0
    d[N-1] = 0.0

    f =  TDMAsolver(a, b, c, d)
    dfdx = np.zeros(N-1)
    for i in range(N-1):
        dfdx[i] = (f[i+1] - f[i]) / dx
		#q[i] = - D[i] * (f[i+1] - f[i]) / dx

    return f, dfdx

# If run as main program
if __name__ == "__main__":

    x_min = 0.0
    x_max = np.pi
    N = 10
    D = np.ones(N-1)
    A = np.ones(N)
    S = np.ones(N)
    # Staggered scheme.
    # Diffusive coefficient.
    xD = np.linspace(x_min, x_max, N-1)
    # Uniform mesh cell size.
    dx = xD[1] - xD[0]
    # Cell centered values.
    x = np.linspace(x_min - dx/2.0, x_max + dx/2.0, N) 

    # Analytic coefficients.
    c0 = 9.
    c1 = 4.
    # Positive diffusion coefficient.
    D = np.sin(xD) + c0
    # Positive ee collisions.
    A = np.sin(x) + c1
    S = -3.0 * np.sin(x) * np.cos(x) - (c0 + c1) * np.cos(x)
    # Analytic solution.
    f_anal = np.cos(x)
    q_anal = - D * (- np.sin(xD)) 

    f, dfdx = IsolatedDiffusionProblem(D, A, S, dx)
    q = - D * dfdx
	#print "min/max f:", min(f), "/", max(f)

    ## Graphics.
    import matplotlib.pyplot as plt
    import matplotlib

    plt.plot(x, f, 'x')
    plt.plot(x, f_anal, label='f')
    plt.legend()
    plt.show()
    plt.plot(xD, q, 'x')
    plt.plot(xD, q_anal, label='q=-Ddfdx')
    plt.legend()
    plt.show()

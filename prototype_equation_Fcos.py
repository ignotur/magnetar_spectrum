from math import *
from scipy.special import legendre
from scipy.optimize import minimize
from scipy import optimize
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rc('font', family='serif')
mpl.rcParams.update({'font.size': 12})
mpl.rcParams.update({'legend.labelspacing':0.25, 'legend.fontsize': 12})
mpl.rcParams.update({'errorbar.capsize': 4})

def fun (alpha, p, N, mu, Pij, d1Pij, d2Pij):

    res = 0.0
    C = alpha[-1]

    res = np.zeros (N+1)

    #print (C)


    for i in range (N):
        #res = res + 
        f = 0.0
        f1 = 0.0
        f2 = 0.0

        for j in range (N):
            f  = f  + alpha[j] * Pij[j,i] 
            f1 = f1 + alpha[j] * d1Pij[j,i]
            f2 = f2 + alpha[j] * d2Pij[j,i]

        res[i] =    (1.0 - mu[i]**2.0)*f2 + p * (p + 1) * f + C * abs(f)**(1.0 + 2.0 / p) 

        if i == 0:
            res[i] = f1


    res[N-2] = f
    res[N-1] = f1 + 2.0

    return res

def f (x, alpha, N):

    res = 0

    for i in range (N):

        #pol   = legendre(i)

        res = res + alpha[i] * cos (i*x) 

    return res

def solve (p, N):

    mu = np.linspace (0.0, 1.0, N)
    alpha = np.zeros (N+1) 
    #alpha[0]  =  2.0/3.0 
    alpha[1]  = 1.0
    alpha[-1] =  0.2


    Pij   = np.zeros ((N,N))
    d1Pij = np.zeros ((N,N))
    d2Pij = np.zeros ((N,N))

    for i in range (0, N):
        #pol   = #legendre(i)
        #pol1d = #pol.deriv(m=1)
        #pol2d = #pol.deriv(m=2)
        for j in range (0, N):

            if i > 0:
                Pij   [i,j] = cos      (i*mu[j])
                d1Pij [i,j] = -i*sin   (i*mu[j])
                d2Pij [i,j] = -i*i*cos (i*mu[j])
            else:
                Pij   [i,j] = 1.0
                d1Pij [i,j] = 0.0
                d2Pij [i,j] = 0.0


    res = optimize.root (fun, alpha, args=(p, N, mu, Pij, d1Pij, d2Pij), method='lm', options={'maxiter':20000})

    print (res)

    return (mu, res.x, res.x[-1])

pl = []
Cl = []

N = 80

mu, alpha, C = solve (0.97, N)

y = []

for i in range (0, N):
    y.append (f (mu[i], alpha, N))

plt.plot (mu, y)
plt.show()





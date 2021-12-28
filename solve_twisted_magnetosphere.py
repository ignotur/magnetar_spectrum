from math import *
import numpy as np
import sys
from scipy import optimize
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rc('font', family='serif')
mpl.rcParams.update({'font.size': 12})
mpl.rcParams.update({'legend.labelspacing':0.25, 'legend.fontsize': 12})
mpl.rcParams.update({'errorbar.capsize': 4})

## Code to solve the non-linear differential equation (6) from Thompson, Lyutikov & Kulkarni (2002)
## The equation is F'' (mu) + p (p + 1) F + C F ^ (1 + 2/p) = 0
## It is solved using the collocation method with cos (nx) basis.
## Boundary conditions are as the following: F'(0) = 0; F (1) = 0 and F'(1) = -2
## C is the parameter which is found together with the solution
## Written by Dr. Andrei P. Igoshev ignotur@gmail.com

## Here we compute the function at the collocation points which are equally spaced 
def fun (alpha, p, N, mu, Pij, d1Pij, d2Pij):

    res = 0.0
    C = alpha[-1]

    res = np.zeros (N+1)

    #print (C)


    for i in range (N):
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

## Here we numerically compute the Jacobian matrix to accelerate the convergance
def dfun_num (alpha, p, N, mu, Pij, d1Pij, d2Pij):

    res = np.zeros ((N+1, N+1))

    for i in range (0, N+1):

        res1 = fun (alpha, p, N, mu, Pij, d1Pij, d2Pij)

        alpha[i] = alpha[i] + 0.01

        res2 = fun (alpha, p, N, mu, Pij, d1Pij, d2Pij)

        alpha[i] = alpha[i] - 0.01

        for j in range (0, N+1):

            res[j,i] = (res2[j] - res1[j]) / 0.01

    return res

## Actual function which initialises and solves the non-linear maxtrix problem using the Levenberg–Marquardt algorithm

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
        for j in range (0, N):

            if i > 0:
                Pij   [i,j] = cos      (i*mu[j])
                d1Pij [i,j] = -i*sin   (i*mu[j])
                d2Pij [i,j] = -i*i*cos (i*mu[j])
            else:
                Pij   [i,j] = 1.0
                d1Pij [i,j] = 0.0
                d2Pij [i,j] = 0.0


    res = optimize.root (fun, alpha, jac=dfun_num, args=(p, N, mu, Pij, d1Pij, d2Pij), method='lm', options={'maxiter':20000})

    print (res)

    return (mu, res.x, res.x[-1])

## To plot function at any mesh
def f (x, alpha, N):
    res = 0
    for i in range (N):
        res = res + alpha[i] * cos (i*x)

    return res



if __name__ == "__main__":

    print ('N = ', sys.argv[1])
    print ('p = ', sys.argv[2])

    try:
        N = int (sys.argv[1])
    except:
        print ('First argument must be integer')
        sys.exit(1)

    try:
        p = float (sys.argv[2])
    except:
        print ('Second argument must be float')
        sys.exit(1)

    mu, alpha, C = solve (p, N)

    y = []

    for i in range (0, N):
        y.append (f (mu[i], alpha, N))

    plt.plot (mu, y)
    plt.xlabel(r'$\mu$')
    plt.ylabel(r'$F(\mu)$')
    plt.savefig ('F_mu.pdf')
    plt.show()

    fl = open ('F_magnetosphere.txt', 'w')

    fl.write (str(N) + '\n')
    fl.write (str(p) + '\n')
    fl.write (str(C) + '\n')

    for i in range (N):
        fl.write (str(mu[i]) + '\t' + str(alpha[i]) + '\n')

    fl.close()








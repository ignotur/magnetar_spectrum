from scipy.integrate import solve_bvp
import sys
import copy
from math import *
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rc('font', family='serif')
mpl.rcParams.update({'font.size': 12})
mpl.rcParams.update({'legend.labelspacing':0.25, 'legend.fontsize': 12})
mpl.rcParams.update({'errorbar.capsize': 4})

def solve (n, p, C, iterat=200, eps_stop = False, eps_lim = 1e-2):

    mu = np.linspace (0, 1.0, n)

    delta_mu = mu[1] - mu[0]

    F = np.sin (pi*mu)
    #F = 1.0 - np.power(mu, 2.0)

    Fnew = copy.copy(F)

    i = 0
    eps = 1
    #for i in range (0, iterat):
    while True:


        if (i > iterat) and (eps_stop == False):
            break
        if (abs(eps) < eps_lim) and eps_stop:
            break

        if eps_stop and i%100 == 0:
            print (i, eps)


        for j in range (1, len(F)-1):

            Fnew [j] = - (1.0-mu[j]**2.0) * (F[j+1] + F[j-1]) / (delta_mu**2.0) / ( p * (p+1) - 2.0 * (1.0-mu[j]**2.0) / (delta_mu**2.0) + C * abs(F[j])**(2.0/p) )

        Fnew[0] = Fnew[1]
        #Fnew[0] = 1.0
        Fnew[n-1] = -2.0 * delta_mu + Fnew[n-2]

        eps = sum((F - Fnew)**2.0)

        F = copy.copy (Fnew)

        i = i + 1

    print (n, p, C, eps)

    return (mu, F)

def f(x,y, C):

    p1 = 0.30
    #C = 0.2

    return np.vstack((y[1], - (C*abs(y[0])**(1.0+2.0/p1) + p1*(p1+1)*y[0] ) / (1.0 - x[0]**2.0)  ))

def bc (ya, yb, C):

    return np.array ([ya[1],yb[0], yb[1]+2])

mu = np.linspace (0.0, 1.0 - 1e-4, 20)

y_a = [np.cos(mu), np.sin(3*mu)]

res_a = solve_bvp(f, bc, mu, y_a, p=[0.2], tol=1e-8)

print (res_a.message)
print (res_a.niter)
print (res_a.p)

y_plot_a = res_a.sol(mu)[0]

plt.plot (mu, y_plot_a)
plt.show()

sys.exit(0)






p = 0.52
C0 = 0.2
n = 20

#mu, F = solve (n, p, 0.2, eps_stop = True, eps_lim = 1e-6)

mu, F = solve (n, p, 0.2, 400)

#print (F[n-1])

print ('-->', F)

plt.plot (mu, F)
plt.show()

sys.exit(0)

while True:

    mu, F = solve (n, p, C0, eps_stop = True, eps_lim = 1e-7)
    C0 = C0 - 100.0 * F[n-1]

    print ('--> ', C0, F[n-1])

    if abs(F[n-1]) < 1e-8:
        break





#mu1, F1 = solve (40, p, 0.01, 200)

#plt.plot (mu1, F1, label='200')

#mu2, F2 = solve (40, p, 0.01, eps_stop = True)

#plt.plot (mu2, F2, label='400')
#plt.legend()
#plt.show()

#C0 = 0.8 

#while True:






sys.exit(0)





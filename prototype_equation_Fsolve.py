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


def f(x,y, C):

    return np.vstack((y[1], - (C*abs(y[0])**(1.0+2.0/p1) + p1*(p1+1.0)*y[0] ) / (1.0 - x[0]**2.0)  ))

def bc (ya, yb, C):

    return np.array ([ya[1],yb[0], yb[1]+2.0])

def solve_for_p1 (pv, f, bc, mu, y_a):

    global p1
    p1 = pv

    res_a = solve_bvp(f, bc, mu, y_a, p=[0.2], tol=1e-8, verbose=0, max_nodes=10000)

    print (res_a.message)
    print (res_a.niter)
    print (res_a.p)

    return [res_a.sol(mu)[0], res_a.p]



mu = np.linspace (0.0, 1.0 , 60)

y_a = [np.cos(mu), np.sin(3*mu)]

y_plot_a01, C01  = solve_for_p1 (0.1, f, bc, mu, y_a)
y_plot_a02, C02  = solve_for_p1 (0.2, f, bc, mu, y_a)
y_plot_a03, C03  = solve_for_p1 (0.3, f, bc, mu, y_a)
y_plot_a04, C04  = solve_for_p1 (0.4, f, bc, mu, y_a)
y_plot_a05, C05  = solve_for_p1 (0.5, f, bc, mu, y_a)
y_plot_a06, C06  = solve_for_p1 (0.6, f, bc, mu, y_a)
y_plot_a07, C07  = solve_for_p1 (0.7, f, bc, mu, y_a)
y_plot_a08, C08  = solve_for_p1 (0.8, f, bc, mu, y_a)
y_plot_a09, C09  = solve_for_p1 (0.9, f, bc, mu, y_a)
y_plot_a10, C10  = solve_for_p1 (1.0, f, bc, mu, y_a)
y_plot_a11, C11  = solve_for_p1 (1.1, f, bc, mu, y_a)


plt.plot (mu, y_plot_a01, 'k-',  label='p = 0.02')
plt.plot (mu, y_plot_a02, 'k--', label='p = 0.2')
plt.plot (mu, y_plot_a03, 'k:',  label='p = 0.3')
plt.plot (mu, y_plot_a05, 'r-',  label='p = 0.5')
plt.plot (mu, y_plot_a10, 'b--', label='p = 1.0')

plt.legend()
plt.show()

pl= [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1]
C = [C01, C02, C03, C04, C05, C06, C07, C08, C09, C10, C11]

plt.scatter (pl, C)
plt.show()








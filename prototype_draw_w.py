from math import *
from scipy.interpolate import interp1d
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rc('font', family='serif')
mpl.rcParams.update({'font.size': 12})
mpl.rcParams.update({'legend.labelspacing':0.25, 'legend.fontsize': 12})
mpl.rcParams.update({'errorbar.capsize': 4})



def f_w (x):
    res = 15.0 * x**2 / pow(pi, 4.0) / (exp(x) - 1)
    return res

def f_n (omega):
    T    = 1e6 ## K
    hbar = 1.05457266e-27 ## erg s
    kB   = 1.380658e-16   ## erg / K
    res = omega * omega / (exp(hbar * omega / kB / T) - 1.0)
    return res

x = np.logspace (-3, 1.8, 100)

fv = []
cf = []

for i in range (0, len(x)):

    fv.append (f_w (x[i]))

cf = np.cumsum (fv) / np.sum(fv)

print (cf)

f = interp1d(cf, x)

plt.plot (x, fv)
plt.xscale('log')
plt.yscale('log')
plt.show()

plt.plot (x, cf)
plt.xscale('log')
#plt.yscale('log')
plt.show()

xx = np.linspace (0.01, 0.99, 100)

plt.plot (cf, x)
plt.plot (xx, f(xx))
plt.show()

T    = 1e6 ## K
hbar = 1.05457266e-27 ## erg s
kB   = 1.380658e-16   ## erg / K

u = np.random.uniform (1e-4, 1.0, size = 5000)

w = []

for i in range (0, len(u)):
    cfv = f(u[i])
    wv =  kB * T / hbar * cfv
#    print (u[i], wv)
    w.append (log10(wv))
#    print (w)
    

omega = np.logspace (15, 19, 100)
omega_log = np.linspace (15, 19, 100)

fn_list = []

for i in range (0, len(omega)):
    f_nv = f_n (omega[i])
    fn_list.append (f_nv)

fn_list = np.asarray(fn_list)
fn_list = fn_list / np.max(fn_list) * 350

plt.hist (w, 50)
plt.plot (omega_log, fn_list)
plt.xlabel(r'$\log_{10} \omega$')
plt.ylabel('Number of photons')
plt.show()



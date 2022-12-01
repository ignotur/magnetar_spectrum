from math import *
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rc('font', family='serif')
mpl.rcParams.update({'font.size': 12})
mpl.rcParams.update({'legend.labelspacing':0.25, 'legend.fontsize': 12})
mpl.rcParams.update({'errorbar.capsize': 4})

f = open ('photon_list.txt', 'r')

E0 = []
E1 = []

for lines in f.readlines ():

    line = lines.split ()

    E0.append (float(line[2]))
    E1.append (float(line[3]))

print ('Minimal initial energy: ', np.min(E0), ' maximum initial energy: ', np.max(E0))
print ('Minimal energy after scattering: ', np.min(E1), ' maximum energy: ', np.max(E1))

bn = np.logspace (-3, 1, 20)

print (bn)

plt.hist (E0, bn, range=[0,10], label='Initial thermal', alpha=0.5)
plt.hist (E1, bn, range=[0,10], label='Scattered', alpha=0.5)
plt.xlabel (r'$E$ (keV)')
plt.ylabel ('Number of photons')
plt.xscale ('log')
plt.yscale ('log')
plt.legend()
plt.savefig('spec.pdf')
plt.show()




from math import *
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rc('font', family='serif')
mpl.rcParams.update({'font.size': 12})
mpl.rcParams.update({'legend.labelspacing':0.25, 'legend.fontsize': 12})
mpl.rcParams.update({'errorbar.capsize': 4})


f = open ('propagation_test.txt', 'r')

x = []; r = []
y = []; theta = []
z = []; phi = []
tau = []

mu = []
omega_fr = []

for lines in f.readlines():

    line = lines.split()

    x.append (float(line[0]))
    y.append (float(line[1]))
    z.append (float(line[2]))

    r.append     (float(line[3]))
    theta.append (float(line[4]))
    phi.append   (float(line[5]))

    tau.append (float(line[6]))

    mu.append (float(line[7]))
    omega_fr.append (float(line[8]))

#plt.plot (r)
#plt.plot (theta)
#plt.plot (phi)
plt.plot (tau)
plt.show()

plt.plot (omega_fr, mu)
plt.xlim([0, 3])
plt.show()




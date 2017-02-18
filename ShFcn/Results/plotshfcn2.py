import numpy as np
import matplotlib.pyplot as plt
from numpy import array

d = np.genfromtxt("xsmall.txt",delimiter=" ")
T = np.genfromtxt("tsmall.txt",delimiter=" ")
Ti = np.genfromtxt("tinterp2.txt",delimiter=" ")
db = np.genfromtxt("xbig.txt",delimiter=" ")
Tb = np.genfromtxt("tbig.txt",delimiter=" ")

# x = array(x_p);
# y = array(T);
axes = plt.gca()
# axes.set_xlim([0,75])
# axes.set_ylim([200,1800])
plt.xlabel('particle location', fontsize=18)
plt.ylabel('temperature (C)', fontsize=16)
# plt.plot((0.001810816, 0.001810816), (200, 1800), 'r-', linewidth=3)

dd = plt.scatter(d, T)
di = plt.scatter(d, Ti, color='green')
db = plt.scatter(db, Tb, color='red')

plt.legend((dd, di, db),
	('Particle level', 'Interpolated', 'Coarse level'))

plt.show()


# axes.set_ylim([0,0003])
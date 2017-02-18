import numpy as np
import matplotlib.pyplot as plt
from numpy import array

d = np.genfromtxt("dist.txt",delimiter=" ")
T = np.genfromtxt("temp.txt",delimiter=" ")

# x = array(x_p);
# y = array(T);
axes = plt.gca()
axes.set_xlim([0,75])
# axes.set_ylim([200,1800])
plt.xlabel('distance from the center node', fontsize=18)
plt.ylabel('temperature (K)', fontsize=16)
# plt.plot((0.001810816, 0.001810816), (200, 1800), 'r-', linewidth=3)

plt.scatter(d, T)
plt.show()


# axes.set_ylim([0,0003])
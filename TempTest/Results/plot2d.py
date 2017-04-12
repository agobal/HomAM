import numpy as np
import matplotlib.pyplot as plt
from numpy import array

x = np.genfromtxt("x.txt",delimiter=" ")
y = np.genfromtxt("T.txt",delimiter=" ")

# x = array(x_p);
# y = array(T);
axes = plt.gca()
axes.set_xlim([0,0.003])
axes.set_ylim([200,1800])
plt.xlabel('distance in x direction (m)', fontsize=18)
plt.ylabel('temperature (K)', fontsize=16)
plt.plot((0.001810816, 0.001810816), (200, 1800), 'r-', linewidth=3)

plt.scatter(x, y)
plt.show()


# axes.set_ylim([0,0003])
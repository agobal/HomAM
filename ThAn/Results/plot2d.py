import numpy as np
import matplotlib.pyplot as plt
from numpy import array

x = np.genfromtxt("x.txt",delimiter=" ")
y = np.genfromtxt("T.txt",delimiter=" ")

# x = array(x_p);
# y = array(T);

plt.scatter(x, y)
plt.show()
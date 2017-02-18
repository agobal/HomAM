import numpy as np
import matplotlib.pyplot as plt
from numpy import array
from matplotlib.legend_handler import HandlerLine2D

xfile = open("runtime1.txt", "r")
rt1 = xfile.read()

yfile = open("runtime2.txt", "r")
rt2 = yfile.read()

zfile = open("numberofparticles.txt", "r")
num = zfile.read()
axes = plt.gca()

rt1 = rt1.split(" ")
rt2 = rt2.split(" ")
num = num.split(" ")

for i, item in enumerate(rt1):
	rt1[i] = float(rt1[i])/1000.0
	rt2[i] = float(rt2[i])/1000.0
	num[i] = int(num[i])

plt.xlabel('number of particles inside packing', fontsize=18)
plt.ylabel('simulation runtime (seconds)', fontsize=16)

pl1 = plt.plot(num, rt1, color='red', label='without grids', linewidth=3)
pl2 = plt.plot(num, rt2, color='blue', label='with grids', linewidth=3)

plt.legend(loc=2)


plt.show()
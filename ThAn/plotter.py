from visual import *

xfile = open("x.txt", "r")
x_p = xfile.read()

yfile = open("y.txt", "r")
y_p = yfile.read()

zfile = open("z.txt", "r")
z_p = zfile.read()

rfile = open("r.txt", "r")
r_p = rfile.read()

x_p = x_p.split(" ")
y_p = y_p.split(" ")
z_p = z_p.split(" ")
r_p = r_p.split(" ")

for idx, val in enumerate(x_p[0:-1]):
	xp = float(x_p[idx])
	yp = float(y_p[idx])
	zp = float(z_p[idx])
	rp = float(r_p[idx])

	sphere(pos=(xp, yp, zp), radius=rp, color=color.cyan)
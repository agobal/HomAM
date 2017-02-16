from visual import *

xfile = open("x.txt", "r")
x_p = xfile.read()

yfile = open("y.txt", "r")
y_p = yfile.read()

zfile = open("z.txt", "r")
z_p = zfile.read()

rfile = open("r.txt", "r")
r_p = rfile.read()

Tfile = open("T.txt", "r")
T = Tfile.read()

x_p = x_p.split(" ")
y_p = y_p.split(" ")
z_p = z_p.split(" ")
r_p = r_p.split(" ")
T = T.split(" ")

for idx, val in enumerate(T):
	try:
		val = float(val)
	except ValueError:
		val = 0
	T[idx] = val

Tmax = 1700
print(max(T))

rmaxidx = len(r_p) - 2
ridx = -1


display(width=800, height=200, center=(0.0012,0.0002,0), background=(1,1,1))

for idx, val in enumerate(x_p[0:-1]):
	# if idx<100:
	xp = float(x_p[idx])
	yp = float(y_p[idx])
	zp = float(z_p[idx])
	ridx = ridx + 1
	if ridx > rmaxidx:
		ridx = 0
	rp = float(r_p[ridx])

	try:
		op = float(T[idx])/float(Tmax)
		if (op > 1):
			print (op)
		# if (op < 0.3):
		# 	op = 0
	except ValueError:
		op = 0
	sphere(pos=(xp, yp, zp), radius=rp, color=color.white)
	if op < 0.25:
		sphere(pos=(xp, yp, zp), radius=rp, color=color.cyan, opacity = op)
	elif op < 0.5:
		sphere(pos=(xp, yp, zp), radius=rp, color=color.green, opacity = op)
	elif op < 0.75:
		sphere(pos=(xp, yp, zp), radius=rp, color=color.yellow, opacity = op)
	else:
		sphere(pos=(xp, yp, zp), radius=rp, color=color.red, opacity = op)
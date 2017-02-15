#include<bits/stdc++.h>
#include "LS.h"
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>

using namespace std;

// Function for setting the initial location of powder particles inside the packing
void PowBed::LaserBeam(int cell, int particle, float speed, int i, float delta_t)
{
	float r_laser = 0.00007;
	Laser_Intensity = ((2*W)/(4.0*atan(1)*r_laser*r_laser))*pow(2.7183, (-(2*(pow((PP.x_p[cell][particle] - i*delta_t*speed), 2) + pow((PP.y_p[cell][particle] - 0.0002), 2)))/(pow(r_laser, 2))));
	Laser_Intensity = Laser_Intensity*exp((-pow((grid_z - PP.z_p[cell][particle]), 2))/(2*pow(2*r_laser, 2)));
}
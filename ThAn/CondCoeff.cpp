#include<bits/stdc++.h>
#include "LS.h"
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>

using namespace std;

// Function for setting the initial location of powder particles inside the packing
void PowBed::CondCoeff(int cell1, int par1, int cell2, int par2, float delta_t)
{
	// This function calculates the heat conduction coefficient between two powder particles based on their location and temperatures
	// Thermal properties
	float r1 = PP.r_p[par1];
	float r2 = PP.r_p[par2];
	float x1 = PP.x_p[cell1][par1];
	float x2 = PP.x_p[cell2][par2];
	float y1 = PP.y_p[cell1][par1];
	float y2 = PP.y_p[cell2][par2];
	float z1 = PP.z_p[cell1][par1];
	float z2 = PP.z_p[cell2][par2];
	float k1 = 22; // Thermal conductivity of particle 1
	float k2 = 22; // Thermal conductivity of particle 2
	float k_air = 0.02; // Thermal conductivity of air
	float alpha1 = 1;
	float alpha2 = 1;
	float k = (2*k1*k2)/(k1 + k2); // Average thermal conductivity of the two surfaces;
	float F = 0.01; // The contact force between the two particles
	float E = 207000000000.0; // Young's modulus of the material
	float nu = 0.33; // Poisson's ratio of the material
	float E_prime = E/(1 - pow(nu, 2));
	float r_eq = (2*r1*r2)/(r1 + r2);
	float m = 0.1; // Average slope of surface roughness
	float sigma = 0.00002; // RMS roughness height
	float r_c, P, A_a, A_r, A_v, cc_dist;
	int s_flag;

	// Mechanical properties
	// if (PP.sintering_flag[cell1][par1] == 0)
	// {
	// 	r_c = pow(((3.0/16.0)*r_eq*F/E_prime) , 0.33); // Contact area radius between two circles
	// 	P = F/(4.0*atan(1)*pow(r_c, 2)); // Compression pressure
	// 	A_a = 4*atan(1)*pow(r_c, 2); // Total contact area
	// 	A_r = A_a*(P*1.41)/(E_prime*m); // Contact point area
	// 	A_v = A_a - A_r; // Void contact area

	// 	cc_dist = r1*cos(asin(r_c/r1)) + r2*cos(asin(r_c/r2));	// Distance of center to center particles
	// 	// Move particle 2 to distance
	// 	PP.x_p[cell2][par2] = x2 + (r1 + r2 - cc_dist)/(r1 + r2)*(x1 - x2);
	// 	PP.y_p[cell2][par2] = y2 + (r1 + r2 - cc_dist)/(r1 + r2)*(y1 - y2);
	// 	PP.z_p[cell2][par2] = z2 + (r1 + r2 - cc_dist)/(r1 + r2)*(z1 - z2);
	// 	PP.sintering_flag[cell2][par2] = 1;
	// }
	// if (PP.sintering_flag[cell1][par1] == 1)
	// {

	// 	cc_dist = sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2));
	// 	r_c = r2*sin(acos(cc_dist/(2*r2)));
	// 	P = F/(4.0*atan(1)*pow(r_c, 2)); // Compression pressure
	// 	A_a = 4*atan(1)*pow(r_c, 2); // Total contact area
	// 	A_r = A_a; // Contact point area
	// 	A_v = A_a - A_r; // Void contact area
	// 	PP.sintering_flag[cell2][par2] = 1;
	// }

	// Material properties or whatever (used for sintering calculations)
	float rho = 2900; // Material density
	float k_n = 700000; // Contact stiffness (mechanical)
	float epsilon = 0.7; // Damping coefficient
	float miu = 0.05; // Friction coefficient
	float D_gDelta_g = 1.8*pow(10, -19); // Diffusion coefficient
	float omega = 1.21*pow(10, -29); // Atomic volume
	float gamma_s = 1.58; // Surface energy
	float psi = (5.0/6.0)*4.0*atan(1); // Dihedral angle
	float k_b = 1.38064852*pow(10, -23); // Boltzmann constant

	float D_b = (D_gDelta_g*omega)/(k_b*PP.T_p[cell1][par1]); // Effective diffusion coefficient
	float v_n = (8*D_b*gamma_s)/(pow(r_c, 4))*(4*r_eq*(1 - cos(psi/2)) + r_c*sin(psi/2));
	// float F_n = (4.0*atan(1)*pow(r_c, 4))/(8*D_b)*v_n - 4.0*atan(1)*gamma_s*(4*r_eq*(1 - cos(psi/2)) + r_c*sin(psi/2));
	float r_c_dot = -r_eq*v_n/r_c; // temporal change in r_c

	r_c = pow(((3.0/16.0)*r_eq*F/E_prime) , 0.33);
	// r_c = r_c + r_c_dot*delta_t;
	A_a = 4*atan(1)*pow(r_c, 2); // Total contact area
	A_r = A_a; // Contact point area
	A_v = A_a - A_r; // Void contact area
	cc_dist = r1*cos(asin(r_c/r1)) + r2*cos(asin(r_c/r2));

	PP.x_p[cell2][par2] = x2 + (r1 + r2 - cc_dist)/(r1 + r2)*(x1 - x2);
	PP.y_p[cell2][par2] = y2 + (r1 + r2 - cc_dist)/(r1 + r2)*(y1 - y2);
	PP.z_p[cell2][par2] = z2 + (r1 + r2 - cc_dist)/(r1 + r2)*(z1 - z2);

	float h_e = (1.55*m*k/sigma)*pow((P*1.41)/(E_prime*m) , 0.94);
	float h_c = h_e + (A_v/(A_a*sigma))*k_air;
	K_F = h_c*A_a; // Overall heat transfer (in room temperature)
}
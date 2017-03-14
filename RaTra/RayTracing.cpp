#include<bits/stdc++.h>
#include "RT.h"
#include <fstream>
#include <math.h>
#include <fstream>

using namespace std;

// Function to read powder bed details from files
void PowderBed::RayTracing()
{
	// Open up files to save x, y, z of ray (for test of correction)
	ofstream xfile("./xx.txt");
	ofstream yfile("./yy.txt");
	ofstream zfile("./zz.txt");
	ofstream pfile("./power.txt");
	ofstream gridfile("./grid.txt");
	ofstream parfile("./par.txt");

	// First we need to find x, y coordinates of the powder bed
	float volume;
	volume = 0;
	for (int i = 1; i <= PP.grid_count; ++i)
	{
		for (int j = 0; j < PP.particle_count; ++j)
		{
			volume = volume + (4.0/3.0)*4.0*atan(1)*pow(PP.r_p[j], 3.0);
		}
	}
	volume = volume/0.57;
	float x_bed = pow(volume, 0.33333);
	float y_bed = x_bed;
	float z_bed = x_bed;

	// Initiating laser ray location randomly
	float x_l, y_l, z_l, theta_l, phi_l, psi_l;

	float A, B, C, delta; // 2nd order equation constants
	float x_2, y_2, z_2, xx, yy, zz, d; // temporary intersection points
	float min_distance;
	float distance;
	int particle_i, particle_gr;
	float vec[3], normal_vec[3], refvec[3];
	float normal_length;
	float cost;
	int prev_grid, prev_par;
	min_distance = 100.0;
	float cosss;
	float lennnn;
	float poww;
	for (int q = 1; q < 10000; ++q)
	{
		x_l = ((double) rand() / (RAND_MAX))*(0.5*x_bed) + 0.25*x_bed;
		y_l = ((double) rand() / (RAND_MAX))*(0.5*y_bed) + 0.25*y_bed;
		z_l = z_bed*2.0;
		theta_l = 2.0*atan(1); // Angle with positive direction of x
		phi_l = 2.0*atan(1); // Angle with positive direction of y;
		psi_l = 4.0*atan(1); // Angle with positive direction of z;
		poww = 1;
		while (poww > 0.01)
		{
			// cout << z_l << " ";
			xfile << x_l << " ";
			yfile << y_l << " ";
			zfile << z_l << " ";
			pfile << poww << " ";
			gridfile << prev_grid << " ";
			parfile << prev_par << " ";
			min_distance = 100.0;
			particle_gr = 0;
			particle_i = 0;
			prev_par = 0;
			prev_grid = 0;
			for (int i = 1; i <= PP.grid_count; ++i)
			{
				for (int j = 0; j < PP.particle_count; ++j)
				{
					if (j == prev_par)
					{
						if (i == prev_grid)
							continue;
						else
							1;
					}
					else
					{
						A = pow(cos(theta_l), 2) + pow(cos(phi_l), 2) + pow(cos(psi_l), 2);
						B = 2*(x_l - PP.x_p[i][j])*cos(theta_l) + 2*(y_l - PP.y_p[i][j])*cos(phi_l) + 2*(z_l - PP.z_p[i][j])*cos(psi_l);
						C = pow(x_l - PP.x_p[i][j], 2) + pow(y_l - PP.y_p[i][j], 2) + pow(z_l - PP.z_p[i][j], 2) - pow(PP.r_p[j], 2);
						delta = pow(B, 2) - 4*A*C;
						lennnn = sqrt(pow((PP.x_p[i][j] - x_l), 2) + pow((PP.y_p[i][j] - y_l), 2) + pow((PP.z_p[i][j] - z_l), 2));
						cosss = cos(theta_l)*(PP.x_p[i][j] - x_l)/lennnn + cos(phi_l)*(PP.y_p[i][j] - y_l)/lennnn + cos(psi_l)*(PP.z_p[i][j] - z_l)/lennnn;
						if ((delta > 0.0) && (cosss >= 0))
						{
							d = (-B + pow(delta, 0.5))/(2*A);
							xx = x_l + d*cos(theta_l);
							yy = y_l + d*cos(phi_l);
							zz = z_l + d*cos(psi_l);
							distance = sqrt(pow(xx - x_l, 2) + pow(yy - y_l, 2) + pow(zz - z_l, 2));
							if (distance < min_distance)
							{
								min_distance = distance;
								x_2 = xx;
								y_2 = yy;
								z_2 = zz;
								particle_gr = i;
								particle_i = j; 
								prev_par = j;
								prev_grid = i;
							}
							d = (-B - pow(delta, 0.5))/(2*A);
							xx = x_l + d*cos(theta_l);
							yy = y_l + d*cos(phi_l);
							zz = z_l + d*cos(psi_l);
							distance = sqrt(pow(xx - x_l, 2) + pow(yy - y_l, 2) + pow(zz - z_l, 2));
							if (distance < min_distance)
							{
								min_distance = distance;
								x_2 = xx;
								y_2 = yy;
								z_2 = zz;
								particle_gr = i;
								particle_i = j; 
								prev_par = j;
								prev_grid = i;
							}
							if (min_distance == 100.0)
								particle_gr = 0;
						}
					}
				}
			}
			if (particle_gr == 0)
				break;
			vec[0] = cos(theta_l);
			vec[1] = cos(phi_l);
			vec[2] = cos(psi_l);

			normal_length = PP.r_p[particle_i];
			// // normal_vec = [(x_2 - PP.x_p[particle_gr][particle_i]), (y_2 - PP.y_p[particle_gr][particle_i]), (z_2 - PP.z_p[particle_gr][particle_i])];
			normal_vec[0] = (x_2 - PP.x_p[particle_gr][particle_i])/normal_length;
			normal_vec[1] = (y_2 - PP.y_p[particle_gr][particle_i])/normal_length;
			normal_vec[2] = (z_2 - PP.z_p[particle_gr][particle_i])/normal_length;

			// Find the angle in between vec and normal vec
			cost = (vec[0]*normal_vec[0] + vec[1]*normal_vec[1] + vec[2]*normal_vec[2])/(sqrt(pow(vec[0], 2) + pow(vec[1], 2) + pow(vec[2], 2))*sqrt(pow(normal_vec[0], 2) + pow(normal_vec[1], 2) + pow(normal_vec[2], 2)));
			float tt = acos(cost)*180.0/(4*atan(1));

			refvec[0] = vec[0] - 2*(vec[0]*normal_vec[0] + vec[1]*normal_vec[1] + vec[2]*normal_vec[2])*normal_vec[0];
			refvec[1] = vec[1] - 2*(vec[0]*normal_vec[0] + vec[1]*normal_vec[1] + vec[2]*normal_vec[2])*normal_vec[1];
			refvec[2] = vec[2] - 2*(vec[0]*normal_vec[0] + vec[1]*normal_vec[1] + vec[2]*normal_vec[2])*normal_vec[2];


			x_l = x_2;
			y_l = y_2;
			z_l = z_2;

			theta_l = acos(refvec[0]);
			phi_l = acos(refvec[1]);
			psi_l = acos(refvec[2]);
			poww = poww*0.7;
			prev_grid = particle_gr;
			prev_par = particle_i;
		}
		xfile << "\n";
		yfile << "\n";
		zfile << "\n";
		pfile << "\n";
		gridfile << "\n";
		parfile << "\n";
		// cout << q << " " << poww << " " << z_l << endl;
	}
	xfile.close();
	yfile.close();
	zfile.close();
	pfile.close();
	gridfile.close();
	parfile.close();
}

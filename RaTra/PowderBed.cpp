#include<bits/stdc++.h>
#include "RT.h"

using namespace std;
 

// Constructor of the arrays needed to perform the simulation on the powder bed structure
PowderBed::PowderBed(int grid_count, int particle_count)
{
	this->particle_count = particle_count;
	this->grid_count = grid_count;

    PP.grid_count = grid_count;
    PP.particle_count = particle_count;
    LR.num_rays = 10000;

    LR.ray_power = new float [LR.num_rays];
    memset(LR.ray_power, 0.0, LR.num_rays*sizeof(float));

    // Create a dynamic array of pointers
    PP.x_p = new float*[grid_count];
    PP.y_p = new float*[grid_count];
    PP.z_p = new float*[grid_count];

    LR.intersection_point = new float*[grid_count];


    // Create a row for every pointer
    for (int i = 1; i <= grid_count; i++)
    {
       // Note : Rows may not be contiguous
       PP.x_p[i] = new float[particle_count];
       PP.y_p[i] = new float[particle_count];
       PP.z_p[i] = new float[particle_count];

       LR.intersection_point[i] = new float[particle_count];

       // Initialize all locations are zero
       memset(PP.x_p[i], 0.0, particle_count*sizeof(float));
       memset(PP.y_p[i], 0.0, particle_count*sizeof(float));
       memset(PP.z_p[i], 0.0, particle_count*sizeof(float));

    }

    // Create the 1D arrays for radii, displacement, force
    PP.r_p = new float [particle_count];
    memset(PP.r_p, 0.0, particle_count*sizeof(float));
}
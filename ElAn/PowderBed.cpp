#include<bits/stdc++.h>
#include "TrussDEM.h"

using namespace std;
 

// Constructor of the arrays needed to perform the simulation on the powder bed structure
PowderBed::PowderBed(int grid_count, int particle_count)
{
	this->particle_count = particle_count;
	this->grid_count = grid_count;

    PP.grid_count = grid_count;
    PP.particle_count = particle_count;

    // Create a dynamic array of pointers
    PP.x_p = new float*[grid_count];
    PP.y_p = new float*[grid_count];
    PP.z_p = new float*[grid_count];
    PP.neighbors = new int**[grid_count];


    // Create a row for every pointer
    for (int i = 1; i <= grid_count; i++)
    {
       // Note : Rows may not be contiguous
       PP.x_p[i] = new float[particle_count];
       PP.y_p[i] = new float[particle_count];
       PP.z_p[i] = new float[particle_count];
       PP.neighbors[i] = new int*[particle_count];

       // Initialize all locations are zero
       memset(PP.x_p[i], 0.0, particle_count*sizeof(float));
       memset(PP.y_p[i], 0.0, particle_count*sizeof(float));
       memset(PP.z_p[i], 0.0, particle_count*sizeof(float));

       for (int j = 0; j < particle_count; ++j)
       {
          PP.neighbors[i][j] = new int[15];
          memset(PP.neighbors[i][j], 0.0, 15*sizeof(int));
        }
    }

    // Create the 1D arrays for radii, displacement, force
    PP.r_p = new float [particle_count];
    memset(PP.r_p, 0.0, particle_count*sizeof(float));


    int par_dof = 6; // DoF for each particle (x, y, z)
    int tot_dof = grid_count*particle_count*par_dof;	// Total displacement and force dof

    PP.u = new float [tot_dof];
    memset(PP.u, 0.0, tot_dof*sizeof(float));

    PP.f = new float [tot_dof];
    memset(PP.f, 0.0, tot_dof*sizeof(float));

    PP.K = new float*[tot_dof];
    PP.M = new float*[tot_dof];
    for (int i = 1; i <= tot_dof; i++)
    {
       // Note : Rows may not be contiguous
       PP.K[i] = new float[tot_dof];
       PP.M[i] = new float[tot_dof];

       // Initialize all locations are zero
       memset(PP.K[i], 0.0, tot_dof*sizeof(float));
       memset(PP.M[i], 0.0, tot_dof*sizeof(float));
    }
}
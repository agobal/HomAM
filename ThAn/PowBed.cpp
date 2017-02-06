#include<bits/stdc++.h>
#include "LS.h"
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>

using namespace std;
 

// Constructor of the arrays needed to perform the simulation on the powder bed structure
PowBed::PowBed(int grid, int par)
{
    this->par = par;
    this->grid = grid;

    PP.grid_count = grid;
    PP.particle_count = par;

    // Create a dynamic array of pointers
    PP.x_p = new float*[grid];
    PP.y_p = new float*[grid];
    PP.z_p = new float*[grid];
    PP.neighbors = new int**[grid];
    PP.sintering_flag = new int*[grid];

    PP.T_p = new float*[grid];
    PP.delta_T_p = new float*[grid];

    // Create a row for every pointer
    for (int i=0; i<grid; i++)
    {
       // Note : Rows may not be contiguous
       PP.x_p[i] = new float[par];
       PP.y_p[i] = new float[par];
       PP.z_p[i] = new float[par];
       PP.neighbors[i] = new int*[par];
       PP.sintering_flag[i] = new int[par];

       PP.T_p[i] = new float[par];
       PP.delta_T_p[i] = new float[par];

       // Initialize all locations are zero
       memset(PP.x_p[i], 0.0, par*sizeof(float));
       memset(PP.y_p[i], 0.0, par*sizeof(float));
       memset(PP.z_p[i], 0.0, par*sizeof(float));
       memset(PP.sintering_flag[i], 0.0, par*sizeof(int));
       memset(PP.T_p[i], 0.0, par*sizeof(float));
       memset(PP.delta_T_p[i], 0.0, par*sizeof(float));

       for (int j = 0; j < par; ++j)
       {
          PP.neighbors[i][j] = new int[15];
          memset(PP.neighbors[i][j], 0.0, 15*sizeof(int));
        }
    }

    PP.r_p = new float [par];
    memset(PP.r_p, 0.0, par*sizeof(float));

}



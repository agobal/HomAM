#include<bits/stdc++.h>
#include "LS.h"
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>

using namespace std;
 

// Constructor of the arrays needed to perform the simulation on the powder bed structure
PowBed::PowBed(int grid, int par, int par_l)
{
    this->par = par;
    this->grid = grid;
    this->par_l = par_l;

    PP.grid_count = grid;
    PP.particle_count = par;

    // For the large packing
    PPl.grid_count = grid;
    PPl.particle_count = par_l;

    // Create a dynamic array of pointers
    PP.x_p = new float*[grid];
    PP.y_p = new float*[grid];
    PP.z_p = new float*[grid];
    PP.neighbors = new int**[grid];
    PP.sintering_flag = new int*[grid];
    // For the large packing
    PPl.x_p = new float*[grid];
    PPl.y_p = new float*[grid];
    PPl.z_p = new float*[grid];
    PPl.neighbors = new int**[grid];
    PPl.sintering_flag = new int*[grid];

    PP.T_p = new float*[grid];
    PP.T_temp = new float*[grid];
    PP.E = new float*[grid];
    // For the large packing
    PPl.T_p = new float*[grid];
    PPl.E = new float*[grid];

    PP.cell_sinter_flag = new int[grid]; // For fast computation (adaptive)
    PP.max_temp = new float[grid];

    // Create a row for every pointer
    for (int i=1; i<=grid; i++)
    {
       // Note : Rows may not be contiguous
       PP.x_p[i] = new float[par];
       PP.y_p[i] = new float[par];
       PP.z_p[i] = new float[par];
       PP.neighbors[i] = new int*[par];
       PP.sintering_flag[i] = new int[par];
       // For the large packing
       PPl.x_p[i] = new float[par_l];
       PPl.y_p[i] = new float[par_l];
       PPl.z_p[i] = new float[par_l];
       PPl.neighbors[i] = new int*[par_l];
       PPl.sintering_flag[i] = new int[par_l];

       PP.T_p[i] = new float[par];
       PP.T_temp[i] = new float[par];
       PP.E[i] = new float[par];
       // For the large packing
       PPl.T_p[i] = new float[par_l];
       PPl.E[i] = new float[par_l];

       // Initialize all locations are zero
       memset(PP.x_p[i], 0.0, par*sizeof(float));
       memset(PP.y_p[i], 0.0, par*sizeof(float));
       memset(PP.z_p[i], 0.0, par*sizeof(float));
       memset(PP.sintering_flag[i], 0.0, par*sizeof(int));
       memset(PP.T_p[i], 0.0, par*sizeof(float));
       memset(PP.T_temp[i], 0.0, par*sizeof(float));
       memset(PP.E[i], 0.0, par*sizeof(float));
       // For the large packing
       memset(PPl.x_p[i], 0.0, par_l*sizeof(float));
       memset(PPl.y_p[i], 0.0, par_l*sizeof(float));
       memset(PPl.z_p[i], 0.0, par_l*sizeof(float));
       memset(PPl.sintering_flag[i], 0.0, par_l*sizeof(int));
       memset(PPl.T_p[i], 0.0, par_l*sizeof(float));
       memset(PPl.E[i], 0.0, par_l*sizeof(float));

       for (int j = 0; j < par; ++j)
       {
          PP.neighbors[i][j] = new int[15];
          memset(PP.neighbors[i][j], 0.0, 15*sizeof(int));
        }
        // Large packing
       for (int j = 0; j < par_l; ++j)
       {
          PPl.neighbors[i][j] = new int[15];
          memset(PPl.neighbors[i][j], 0.0, 15*sizeof(int));
        }
    }

    PP.r_p = new float [par];
    memset(PP.r_p, 0.0, par*sizeof(float));
    // Large packing    
    PPl.r_p = new float [par_l];
    memset(PPl.r_p, 0.0, par_l*sizeof(float));

}



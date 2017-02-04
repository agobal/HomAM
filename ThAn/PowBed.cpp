#include<bits/stdc++.h>
#include "LS.h"
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>

using namespace std;
 

 
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
    PP.sintering_flag = new int*[grid];

    // Create a row for every pointer
    for (int i=0; i<grid; i++)
    {
       // Note : Rows may not be contiguous
       PP.x_p[i] = new float[par];
       PP.y_p[i] = new float[par];
       PP.z_p[i] = new float[par];
       PP.sintering_flag[i] = new int[par];

       // Initialize all locations are zero
       memset(PP.x_p[i], 0.0, par*sizeof(float));
       memset(PP.y_p[i], 0.0, par*sizeof(float));
       memset(PP.z_p[i], 0.0, par*sizeof(float));
       memset(PP.sintering_flag[i], 0.0, par*sizeof(int));
    }

    PP.r_p = new float [par];
    memset(PP.r_p, 0.0, par*sizeof(float));

}
 
// Utility method to print adjacency matrix
void PowBed::print()
{
   for (int u=0; u<PP.particle_count; u++)
   {
      // for (int v=0; v<PP.particle_count; v++)
         cout << PP.r_p[u] << " ";
      cout << endl;
   }
}

// Set initial locations
void PowBed::SetInitials(int grid, int par)
{
  num_grid_x = 15;
  num_grid_y = 15;
  num_grid_z = 15;

  PP.avgrd = 0.0125;
  PP.stddev = 0.005;

  /* First create the size distribution using the boost library */
  boost::mt19937 rng;

  boost::normal_distribution<> nd(PP.avgrd, PP.stddev);

  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_nor(rng, nd);

  float volume; // Total volume of the grid cell occupied by the powder particles
  volume = 0;
  float r;        // Intermediate variable for particle radius

  /* Adding particle volumes to find the volume of the powder bed */
  for (int i = 0; i <= par; ++i)
  {
    r = var_nor();
    volume = volume + (4.0/3.0)*4.0*atan(1.0)*pow(r, 3.0);
    // If the volume is cell than the maximum packing volume, add more
    PP.r_p[i] = r;
  }

  // Small grid sizes are determined as:
  grid_x = pow(volume, 0.3);
  grid_y = grid_x;
  grid_z = grid_y;

  /* Randomly putting particles inside the powder bed */

  int cell_x_num[grid], cell_y_num[grid], cell_z_num[grid]; // Number of the cell grid in x, y and z directions

  for (int i = 0; i < grid; ++i)
  {
    // Assigning cell grid numbering (starts from 1 unlike the total cell number which starts from 0)
    cell_x_num[i] = (i % num_grid_x) + 1;
    cell_y_num[i] = (i % (num_grid_x*num_grid_y))/num_grid_x + 1;
    cell_z_num[i] = floor(i/(num_grid_x*num_grid_y) ) + 1;

    for (int j = 0; j < par; ++j)
    {
      // Assigning initial random locations to particles
      PP.x_p[i][j] = ((double) rand() / (RAND_MAX))*(grid_x - 2.0*PP.r_p[i]) + PP.r_p[i] + float((cell_x_num[i] - 1))*grid_x;
      PP.y_p[i][j] = ((double) rand() / (RAND_MAX))*(grid_y - 2.0*PP.r_p[i]) + PP.r_p[i] + float((cell_y_num[i] - 1))*grid_y;
      PP.z_p[i][j] = ((double) rand() / (RAND_MAX))*(grid_z - 2.0*PP.r_p[i]) + PP.r_p[i] + float((cell_z_num[i] - 1))*grid_z;
      PP.sintering_flag[i][j] = 0;
    }
    cout << i << endl;
  }
}

void PowBed::PowDel()
{
   for (int i=0; i<PP.grid_count; i++)
    {
      delete[] PP.x_p[i];
      delete[] PP.y_p[i];
      delete[] PP.z_p[i];
    }
   delete[] PP.x_p;
   delete[] PP.y_p;
   delete[] PP.z_p;
   delete[] PP.r_p;
 }
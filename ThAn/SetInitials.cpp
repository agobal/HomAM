#include<bits/stdc++.h>
#include "LS.h"
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>

using namespace std;

// Function for setting the initial location of powder particles inside the packing
void PowBed::SetInitials(int grid, int par)
{
  num_grid_x = 4;
  num_grid_y = 4;
  num_grid_z = 4;


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

  // Rearranging particles
  float q, x_particle_middle, y_particle_middle, z_particle_middle; // Middle variables for ease of calculations
  int counter, neighbor_part, neighbor_cell;    // loop counter and variables for determining particle neighbors

  // Relocate particles to reduce overlaps
  for (int c = 0; c < grid; ++c)
  {
    cout << cell_x_num[c] << " " << cell_y_num[c] << " " << cell_z_num[c] << endl << endl;
    cout << c << endl;
    // Filling out the neighboring particle numbers for relocation purposes
    // Neighboring list consists of 4 cell groups for the cell itself plus neighbors in x, y, z directions behind it
    int neighbor_particles[4*par];
    counter = 0;
    // Initiating the neighboring particles array
    for (int ct = 0; ct < 4*par; ++ct)
      neighbor_particles[ct] = 0;
    // Adding the particles inside cell to the neighboring list
    for (int c1 = 0; c1 < par; ++c1)
    {
      neighbor_particles[counter] = c1;
      counter = counter + 1;
    }
    // Adding the particles in the cell behind to neighboring list (x)
    if ((num_grid_x != 1) && (cell_x_num[c] != 1))
    {
      for (int c1 = 0; c1 < par; ++c1)
      {
        neighbor_particles[counter] = c1 + (c)*1000; // +1 because the total cell number starts from 0 but 0*1000=0 so the difference wouldn't be obvious therefore we add 1 then subtract it later on
        counter = counter + 1;
      }
    }
    // Adding the particles in the cell behind to neighboring list (y)
    if ((num_grid_y != 1) && (cell_y_num[c] != 1))
    {
      for (int c1 = 0; c1 < par; ++c1)
      {
        neighbor_particles[counter] = c1 + (c)*1000;
        counter = counter + 1;
      }
    }
    // Adding the particles in the cell behind to neighboring list (z)
    if ((num_grid_z != 1) && (cell_z_num[c] != 1))
    {
      for (int c1 = 0; c1 < par; ++c1)
      {
        neighbor_particles[counter] = c1 + (c)*1000;
        counter = counter + 1;
      }
    }
/*    for (int ct = 0; ct < 4*ppcount; ++ct)
      cout << neighbor_particles[ct] << " ";
    cout << endl;*/   //Uncomment for testing particle locations

    // Performing the relocation of particles to reach a stable position
    for (int k = 0; k < 600; ++k)
    {
      for (int i = 0; i < par; ++i)
      {
        q = 0.0;
        x_particle_middle = 0.0;
        y_particle_middle = 0.0;
        z_particle_middle = 0.0;
        for (int j = 0; j < 4*par; ++j)
        {
          // For particles in the same cell
          if (j < par)
          {
            if (j != i)
            {
              if ((PP.r_p[i] + PP.r_p[j]) > sqrt(pow(PP.x_p[c][i] - PP.x_p[c][j], 2.0) + pow(PP.y_p[c][i] - PP.y_p[c][j], 2.0) + pow(PP.z_p[c][i] - PP.z_p[c][j], 2.0)))
              {

                q = q + 1;
                x_particle_middle = x_particle_middle + PP.x_p[c][j] + (PP.x_p[c][i] - PP.x_p[c][j])*(PP.r_p[i] + PP.r_p[j])/sqrt(pow((PP.x_p[c][i] - PP.x_p[c][j]),2.0) + pow((PP.y_p[c][i] - PP.y_p[c][j]),2.0) + pow((PP.z_p[c][i] - PP.z_p[c][j]),2.0));
                y_particle_middle = y_particle_middle + PP.y_p[c][j] + (PP.y_p[c][i] - PP.y_p[c][j])*(PP.r_p[i] + PP.r_p[j])/sqrt(pow((PP.x_p[c][i] - PP.x_p[c][j]),2.0) + pow((PP.y_p[c][i] - PP.y_p[c][j]),2.0) + pow((PP.z_p[c][i] - PP.z_p[c][j]),2.0));
                z_particle_middle = z_particle_middle + PP.z_p[c][j] + (PP.z_p[c][i] - PP.z_p[c][j])*(PP.r_p[i] + PP.r_p[j])/sqrt(pow((PP.x_p[c][i] - PP.x_p[c][j]),2.0) + pow((PP.y_p[c][i] - PP.y_p[c][j]),2.0) + pow((PP.z_p[c][i] - PP.z_p[c][j]),2.0));
              }
            }
          }
          // For particles in the adjacent cells
          else
          {
            neighbor_cell = (neighbor_particles[j]/1000); // Because the total cell count starts from 0 and xyz count starts from 1
            // if ((k == 0) && (i == 0) && (j == par))
            //   cout << neighbor_cell << endl;
            if (neighbor_cell != 0)
              neighbor_part = (neighbor_particles[j] % (neighbor_cell*1000));
            // cout << neighbor_cell << " " << neighbor_part << endl;
            if ((neighbor_cell != 0) && (neighbor_part != 0))
            {
              if ((PP.r_p[i] + PP.r_p[neighbor_part]) > sqrt(pow(PP.x_p[c][i] - PP.x_p[neighbor_cell][neighbor_part], 2.0) + pow(PP.y_p[c][i] - PP.y_p[neighbor_cell][neighbor_part], 2.0) + pow(PP.z_p[c][i] - PP.z_p[neighbor_cell][neighbor_part], 2.0)))
              {
                q = q + 1;
                x_particle_middle = x_particle_middle + PP.x_p[neighbor_cell][neighbor_part] + (PP.x_p[c][i] - PP.x_p[neighbor_cell][neighbor_part])*(PP.r_p[i] + PP.r_p[neighbor_part])/sqrt(pow((PP.x_p[c][i] - PP.x_p[neighbor_cell][neighbor_part]),2.0) + pow((PP.y_p[c][i] - PP.y_p[neighbor_cell][neighbor_part]),2.0) + pow((PP.z_p[c][i] - PP.z_p[neighbor_cell][neighbor_part]),2.0));
                y_particle_middle = y_particle_middle + PP.y_p[neighbor_cell][neighbor_part] + (PP.y_p[c][i] - PP.y_p[neighbor_cell][neighbor_part])*(PP.r_p[i] + PP.r_p[neighbor_part])/sqrt(pow((PP.x_p[c][i] - PP.x_p[neighbor_cell][neighbor_part]),2.0) + pow((PP.y_p[c][i] - PP.y_p[neighbor_cell][neighbor_part]),2.0) + pow((PP.z_p[c][i] - PP.z_p[neighbor_cell][neighbor_part]),2.0));
                z_particle_middle = z_particle_middle + PP.z_p[neighbor_cell][neighbor_part] + (PP.z_p[c][i] - PP.z_p[neighbor_cell][neighbor_part])*(PP.r_p[i] + PP.r_p[neighbor_part])/sqrt(pow((PP.x_p[c][i] - PP.x_p[neighbor_cell][neighbor_part]),2.0) + pow((PP.y_p[c][i] - PP.y_p[neighbor_cell][neighbor_part]),2.0) + pow((PP.z_p[c][i] - PP.z_p[neighbor_cell][neighbor_part]),2.0));
              }
            }
          }
        }
        if (q >= 1)
        {
          PP.x_p[c][i] = x_particle_middle/q;
          PP.y_p[c][i] = y_particle_middle/q;
          PP.z_p[c][i] = z_particle_middle/q;

          // Particles can't go over the overall boundaries of the PB
          if ((cell_x_num[c] == num_grid_x) && (PP.x_p[c][i] >= (grid_x*cell_x_num[c] - PP.r_p[i])))
            PP.x_p[c][i] = (grid_x - PP.r_p[i]);
          if ((cell_y_num[c] == num_grid_y) && (PP.y_p[c][i] >= (grid_y*cell_y_num[c] - PP.r_p[i])))
            PP.y_p[c][i] = (grid_y - PP.r_p[i]);
          if ((cell_z_num[c] == num_grid_z) && (PP.z_p[c][i] >= (grid_z*cell_z_num[c] - PP.r_p[i])))
           PP.z_p[c][i] = (grid_z - PP.r_p[i]);
          // Particles can't go below a certain amount of the previous cell
          if (PP.x_p[c][i] <= (0.9*grid_x*(cell_x_num[c] - 1.0) + PP.r_p[i]))
            PP.x_p[c][i] = (0.9*grid_x*(cell_x_num[c] - 1.0) + PP.r_p[i]);
          if (PP.y_p[c][i] <= (0.9*grid_y*(cell_y_num[c] - 1.0) + PP.r_p[i]))
            PP.y_p[c][i] = (0.9*grid_y*(cell_y_num[c] - 1.0) + PP.r_p[i]);
          if (PP.z_p[c][i] <= (0.9*grid_z*(cell_z_num[c] - 1.0) + PP.r_p[i]))
            PP.z_p[c][i] = (0.9*grid_z*(cell_z_num[c] - 1.0) + PP.r_p[i]);
          // Particles can't go below the general size of the PB
          if ((cell_x_num[c] == 1) && (PP.x_p[c][i] <= (grid_x*(cell_x_num[c] - 1) + PP.r_p[i])))
            PP.x_p[c][i] = PP.r_p[i];
          if ((cell_y_num[c] == 1) && (PP.y_p[c][i] <= (grid_y*(cell_y_num[c] - 1) + PP.r_p[i])))
            PP.y_p[c][i] = PP.r_p[i];
          if ((cell_z_num[c] == 1) && (PP.z_p[c][i] <= (grid_z*(cell_z_num[c] - 1) + PP.r_p[i])))
            PP.z_p[c][i] = PP.r_p[i];
        }
      }
    }
    // Find neighbors of particles
    int q2;
    for (int i = 0; i < par; ++i)
    {
      q2 = 0;
      for (int j = 0; j < 4*par; ++j)
      {
        /// For particles in the same cell
        if (j < par)
        {
          if (j != i)
          {
            if (fabs(PP.r_p[i] + PP.r_p[j] - pow(pow(PP.x_p[c][i] - PP.x_p[c][j], 2.0) + pow(PP.y_p[c][i] - PP.y_p[c][j], 2.0) + pow(PP.z_p[c][i] - PP.z_p[c][j], 2.0), 0.5)) < 0.000001)
            {
              PP.neighbors[c][i][q2] = j;
              q2 = q2 + 1;
            }
          }
        }
        // For particles in the adjacent cells
        else
        {
          neighbor_cell = (neighbor_particles[j]/1000); // Because the total cell count starts from 0 and xyz count starts from 1
          if (neighbor_cell != 0)
            neighbor_part = (neighbor_particles[j] % (neighbor_cell*1000));
          if ((neighbor_cell != 0) && (neighbor_part != 0))
          {
            if (fabs(PP.r_p[i] + PP.r_p[neighbor_part] - sqrt(pow(PP.x_p[c][i] - PP.x_p[neighbor_cell][neighbor_part], 2.0) + pow(PP.y_p[c][i] - PP.y_p[neighbor_cell][neighbor_part], 2.0) + pow(PP.z_p[c][i] - PP.z_p[neighbor_cell][neighbor_part], 2.0))) < 0.000001)
            {
              PP.neighbors[c][i][q2] = neighbor_particles[j];
              q2 = q2 + 1;
            }
          }
        }
      }
    }
  }
}
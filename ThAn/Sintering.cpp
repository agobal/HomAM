#include<bits/stdc++.h>
#include "LS.h"
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>

using namespace std;

// Function for setting the initial location of powder particles inside the packing
void PowBed::Sintering(float power, float speed)
{
  W = power;
  dt = (grid_y*num_grid_y/speed)/1000.0;

  int Q; // Total heat delivered to each particle at each point in time
  int neigh; // intermediate variable for neighbor value

  // Initial temperature of powder particles
  for (int cell = 0; cell < grid; ++cell)
  {
    for (int particle = 0; particle < par; ++particle)
    {
      PP.T_p[cell][particle] = 100;
    }
  }

  // Go through the loop for performing the laser sintering (in the y direction)
  for (int i = 0; i < 1000; ++i)
  {
    for (int cell = 0; cell < grid; ++cell)
    {
      for (int particle = 0; particle < par; ++particle)
      {
        Q = 0;
        // First step in to calculate heat transfer to and from neighbors from and to the designated particle
        for (int j = 0; j < 15; ++j)
        {
          if (PP.neighbors[cell][particle][j] != 0)
          {
            if (PP.neighbors[cell][particle][j] < 1000)
            {
              neigh = PP.neighbors[cell][particle][j];
              // function to calculate the conduction coefficient between powder particles
              CondCoeff(cell, particle, cell, neigh, dt);
              // Calculate the heat transfer between neighbors
              Q = Q + K_F*(PP.T_p[cell][particle] - PP.T_p[cell][neigh]);
            }
            else
            {
              int cell2 = (PP.neighbors[cell][particle][j]/1000);
              if (cell2 != 0)
                neigh = (PP.neighbors[cell][particle][j] % (cell2*1000));
              CondCoeff(cell, particle, cell2, neigh, dt);
              Q = Q + K_F*(PP.T_p[cell][particle] - PP.T_p[cell2][neigh]);
            }
          }
          // Heat input from the laser beam
          LaserBeam(cell, particle, speed, i, dt);

          float S = 4.0*atan(1)*PP.r_p[particle]*PP.r_p[particle];    // Particle surface absorbing the laser powder
          float rho = 7800;
          float K_ab = 0.3;
          float solid_heat_capacity = 477;
          PP.E[cell][particle] = PP.E[cell][particle] + (Q + K_ab*S*Laser_Intensity*dt/(rho*(4.0/3.0)*4.0*atan(1)*pow(PP.r_p[particle], 3))); //particle energy increase by laser
          PP.T_temp[cell][particle] = PP.E[cell][particle]/solid_heat_capacity; // Particle temperature change
        }
      }
    }
    // Populating the temperature array
    for (int cell = 0; cell < grid; ++cell)
    {
      for (int particle = 0; particle < par; ++particle)
      {
        PP.T_p[cell][particle] = PP.T_temp[cell][particle];
      }
    }
  }
}
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
  dt = (grid_x*num_grid_x/speed)/1000.0;
  cout << dt << endl;
  float Q; // Total heat delivered to each particle at each point in time
  int neigh; // intermediate variable for neighbor value

  float rho = 7800;
  float K_ab = 0.3;
  float solid_heat_capacity = 477;

  // Initial temperature of powder particles
  for (int cell = 1; cell <= grid; ++cell)
  {
    for (int particle = 0; particle < par; ++particle)
    {
      PP.T_p[cell][particle] = 10;
      PP.E[cell][particle] = PP.T_p[cell][particle]*(rho*(4.0/3.0)*4.0*atan(1)*pow(PP.r_p[particle], 3))*solid_heat_capacity;
    }
  }

  // Go through the loop for performing the laser sintering (in the y direction)
  for (int i = 0; i < 10000; ++i)
  {
    for (int cell = 1; cell <= grid; ++cell)
    {
      // Boundary condition
      if (cell == 1)
      {
        for (int p = 0; p < par; ++p)
        {
          PP.T_p[cell][p] = 100;
          PP.E[cell][p] = PP.T_p[cell][p]*(rho*(4.0/3.0)*4.0*atan(1)*pow(PP.r_p[p], 3))*solid_heat_capacity;
        }
      }
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
              Q = Q - K_F*(PP.T_p[cell][particle] - PP.T_p[cell][neigh]);
              // if ((cell == 6) && (particle == 10))
              //   cout << K_F << " " << cell << " " << particle << " " << neigh << " " << PP.T_p[cell][particle] << " " << PP.T_p[cell][neigh] << " " << Q << " " << PP.E[cell][particle] << " " << endl;
            }
            else
            {
              int cell2 = (PP.neighbors[cell][particle][j]/1000);
              if (cell2 != 0)
                neigh = (PP.neighbors[cell][particle][j] % (cell2*1000));
              CondCoeff(cell, particle, cell2, neigh, dt);
              Q = Q - K_F*(PP.T_p[cell][particle] - PP.T_p[cell2][neigh]);
            }
          }
        }
        // Heat input from the laser beam
        // LaserBeam(cell, particle, speed, i, dt);
        // + K_ab*S*Laser_Intensity*dt
        float S = 4.0*atan(1)*PP.r_p[particle]*PP.r_p[particle];    // Particle surface absorbing the laser powder
        PP.E[cell][particle] = PP.E[cell][particle] + (Q*dt); //particle energy increase by laser
        PP.T_p[cell][particle] = PP.E[cell][particle]/(solid_heat_capacity*rho*(4.0/3.0)*4.0*atan(1)*pow(PP.r_p[particle], 3)); // Particle temperature change
        // if ((cell == 6) && (particle == 10))
        //   cout << Q << " " << PP.E[cell][particle] << " " << PP.T_p[cell][particle] << " " << (solid_heat_capacity*rho*(4.0/3.0)*4.0*atan(1)*pow(PP.r_p[particle], 3)) << endl;
      }
    }
    // Populating the temperature array
    // for (int cell = 1; cell <= grid; ++cell)
    // {
    //   for (int particle = 0; particle < par; ++particle)
    //   {
    //     PP.T_p[cell][particle] = PP.T_temp[cell][particle];
    //   }
    // }
  }
}
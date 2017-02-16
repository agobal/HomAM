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
  dt = (grid_x*num_grid_x/speed)/10000.0;
  cout << dt << endl;
  float Q; // Total heat delivered to each particle at each point in time
  int neigh; // intermediate variable for neighbor value

  float rho = 7800;
  float K_ab = 0.3;
  float solid_heat_capacity = 500;

  // Initial temperature of powder particles
  for (int cell = 1; cell <= grid; ++cell)
  {
    cout << cell << endl;
    for (int particle = 0; particle < par; ++particle)
    {
      PP.T_p[cell][particle] = 0;
      PP.E[cell][particle] = 0;
    }
    for (int particles = 0; particles < par_l; ++particles)
    {
      if ((cell == 1) && (particles == 10))
        PP.T_p[cell][PPl.sintering_flag[cell][particles]] = 1;
    }
    for (int particle = 0; particle < par_l; ++particle)
    {
      PP.E[cell][PPl.sintering_flag[cell][particle]] = PP.T_p[cell][PPl.sintering_flag[cell][particle]]*(rho*(4.0/3.0)*4.0*atan(1)*pow(PP.r_p[particle], 3))*solid_heat_capacity;
    }
  }

  float xc, yc, zc, dc;
  int cell_x_num[grid], cell_y_num[grid], cell_z_num[grid]; // Number of the cell grid in x, y and z directions
  // Go through the loop for performing the laser sintering (in the y direction)
  for (int i = 0; i < 20000; ++i)
  {
    // cout << i << endl;
    for (int cell = 1; cell <= grid; ++cell)
    {
      PP.max_temp[cell] = 0;
      // Assigning cell grid numbering (starts from 1 unlike the total cell number which starts from 0)
      cell_x_num[cell] = (cell % num_grid_x);
      if (cell_x_num[cell] == 0)
        cell_x_num[cell] = num_grid_x;
      cell_y_num[cell] = (cell % (num_grid_y*num_grid_x));
      if (cell_y_num[cell] == 0)
        cell_y_num[cell] = num_grid_y*num_grid_x;
      cell_y_num[cell] = (cell_y_num[cell] - 1)/num_grid_x + 1;
      cell_z_num[cell] = ((i - 1) / (num_grid_x*num_grid_y)) + 1;
      xc = (cell_x_num[cell]*grid_x + (cell_x_num[cell] - 1)*grid_x)/2;
      yc = (cell_y_num[cell]*grid_y + (cell_y_num[cell] - 1)*grid_y)/2;
      zc = (cell_z_num[cell]*grid_z + (cell_z_num[cell] - 1)*grid_z)/2;
      dc = sqrt(pow(xc - i*dt*speed, 2) + pow(yc - 0.0005, 2));

      for (int particle = 0; particle < par; ++particle)
      {
        Q = 0;
        if (sqrt(pow(PP.x_p[cell][particle] - PP.x_p[cell][PPl.sintering_flag[cell][10]], 2) + pow(PP.y_p[cell][particle] - PP.y_p[cell][PPl.sintering_flag[cell][10]], 2) + pow(PP.z_p[cell][particle] - PP.z_p[cell][PPl.sintering_flag[cell][10]], 2)) > 0.00009)
          PP.T_p[cell][particle] = 0;
        if (particle < par_l)
        {
          if ((cell == 1) && (particle == 10))
          {
            PP.T_p[cell][PPl.sintering_flag[cell][particle]] = 1;
          }
          else
          {
            PP.T_p[cell][PPl.sintering_flag[cell][particle]] = 0;
          }
          PP.E[cell][PPl.sintering_flag[cell][particle]] = PP.T_p[cell][PPl.sintering_flag[cell][particle]]*(rho*(4.0/3.0)*4.0*atan(1)*pow(PP.r_p[PPl.sintering_flag[cell][particle]], 3))*solid_heat_capacity;
        }
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
              Q = Q - 0.001*(PP.T_p[cell][particle] - PP.T_p[cell][neigh]);
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

        float S = 4.0*atan(1)*PP.r_p[particle]*PP.r_p[particle];    // Particle surface absorbing the laser powder
        PP.E[cell][particle] = PP.E[cell][particle] + (Q*dt); //particle energy increase by laser

        // Temperature and melting calculations:
        PP.T_p[cell][particle] = PP.E[cell][particle]/(solid_heat_capacity*rho*(4.0/3.0)*4.0*atan(1)*pow(PP.r_p[particle], 3)); // Particle temperature change

        if (PP.T_p[cell][particle] > PP.max_temp[cell])
          PP.max_temp[cell] = PP.T_p[cell][particle];

        // Assign temperatures to the larger particles
        if (particle < par_l)
        {
          PPl.T_p[cell][particle] = PP.T_p[cell][PPl.sintering_flag[cell][particle]];
          PPl.E[cell][particle] = 15*PPl.T_p[cell][particle]*(rho*(4.0/3.0)*4.0*atan(1)*pow(PP.r_p[particle], 3))*solid_heat_capacity;
        }
      }
      // for (int particle = 0; particle < par_l; ++particle)
      // {
      //   Q = 0;
      //   // First step in to calculate heat transfer to and from neighbors from and to the designated particle
      //   for (int j = 0; j < 15; ++j)
      //   {
      //     if (PPl.neighbors[cell][particle][j] != 0)
      //     {
      //       if (PPl.neighbors[cell][particle][j] < 1000)
      //       {
      //         neigh = PPl.neighbors[cell][particle][j];
      //         // function to calculate the conduction coefficient between powder particles
      //         // CondCoeff(cell, particle, cell, neigh, dt);
      //         // Calculate the heat transfer between neighbors
      //         Q = Q - 10*K_F*(PPl.T_p[cell][particle] - PPl.T_p[cell][neigh]);
      //         // if ((cell == 6) && (particle == 10))
      //         //   cout << K_F << " " << cell << " " << particle << " " << neigh << " " << PP.T_p[cell][particle] << " " << PP.T_p[cell][neigh] << " " << Q << " " << PP.E[cell][particle] << " " << endl;
      //       }
      //       else
      //       {
      //         int cell2 = (PPl.neighbors[cell][particle][j]/1000);
      //         if (cell2 != 0)
      //           neigh = (PPl.neighbors[cell][particle][j] % (cell2*1000));
      //         // CondCoeff(cell, particle, cell2, neigh, dt);
      //         Q = Q - 10*K_F*(PPl.T_p[cell][particle] - PPl.T_p[cell2][neigh]);
      //       }
      //     }
      //   }
      //   // Heat input from the laser beam
      //   PPl.E[cell][particle] = PPl.E[cell][particle] + ((Q)*dt); //particle energy increase by laser
      //   PPl.T_p[cell][particle] = PPl.E[cell][particle]/(0.5*solid_heat_capacity*rho*(4.0/3.0)*4.0*atan(1)*pow(PPl.r_p[particle], 3)); // Particle temperature change
      // }
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
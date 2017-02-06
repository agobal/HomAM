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
              // function to calculated the conduction coefficient between powder particles
              // Calculate the heat transfer between neighbors
              // Calculate the  heat input from the laser beam
              // Calculate displacement of powder particles from the sintering degree (dozaria...)
            }
          }
        }
      }
    }
  }
}
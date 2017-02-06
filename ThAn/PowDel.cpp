#include<bits/stdc++.h>
#include "LS.h"
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>

using namespace std;

// Function for clearing up the heap memory after dealing with the problem
void PowBed::PowDel()
{
   for (int i=0; i<PP.grid_count; i++)
    {
      delete[] PP.x_p[i];
      delete[] PP.y_p[i];
      delete[] PP.z_p[i];
      delete[] PP.T_p[i];
      delete[] PP.delta_T_p[i];
    }
   delete[] PP.x_p;
   delete[] PP.y_p;
   delete[] PP.z_p;
   delete[] PP.r_p;
   delete[] PP.T_p;
   delete[] PP.delta_T_p;
 }
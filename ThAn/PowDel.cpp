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
   for (int i=1; i<=PP.grid_count; i++)
    {
      delete[] PP.x_p[i];
      delete[] PP.y_p[i];
      delete[] PP.z_p[i];
      delete[] PP.T_p[i];
      delete[] PP.T_temp[i];
      delete[] PP.E[i];
      delete[] PP.sintering_flag[i];

      //LArge particles
      delete[] PPl.x_p[i];
      delete[] PPl.y_p[i];
      delete[] PPl.z_p[i];
      delete[] PPl.T_p[i];
      delete[] PPl.E[i];
    }
   delete[] PP.x_p;
   delete[] PP.y_p;
   delete[] PP.z_p;
   delete[] PP.r_p;
   delete[] PP.T_p;
   delete[] PP.T_temp;
   delete[] PP.E;
   delete[] PP.sintering_flag;
   delete[] PP.cell_sinter_flag;

   // Large particles
   delete[] PPl.x_p;
   delete[] PPl.y_p;
   delete[] PPl.z_p;
   delete[] PPl.r_p;
   delete[] PPl.T_p;
   delete[] PPl.E;
 }
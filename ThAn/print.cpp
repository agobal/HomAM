#include<bits/stdc++.h>
#include "LS.h"
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>

using namespace std;

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

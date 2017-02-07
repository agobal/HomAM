#include<bits/stdc++.h>
#include "LS.h"
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <iostream>
#include <fstream>

using namespace std;

// Utility method to print adjacency matrix
void PowBed::print()
{
	ofstream xfile("x.txt");
	for (int i = 0; i < PP.grid_count; ++i)
	{
		for (int j = 0; j < PP.particle_count; ++j)
		{
			xfile << PP.x_p[i][j];
			xfile << " ";
		}
		xfile << "\n";
	}
	xfile.close();
   // for (int u=0; u<PP.particle_count; u++)
   // {
   //   cout << PP.r_p[u] << " ";
   // }
}

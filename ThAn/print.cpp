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
	ofstream yfile("y.txt");
	ofstream zfile("z.txt");
	ofstream rfile("r.txt");

	ofstream Tfile("T.txt");

	ofstream Neighfile("Neighbors.txt");

	for (int i = 0; i < PP.particle_count; ++i)
	{
		rfile << PP.r_p[i];
		rfile << " ";
	}
	for (int i = 0; i < PP.grid_count; ++i)
	{
		for (int j = 0; j < PP.particle_count; ++j)
		{
			xfile << PP.x_p[i][j];
			xfile << " ";
			yfile << PP.y_p[i][j];
			yfile << " ";
			zfile << PP.z_p[i][j];
			zfile << " ";
			Tfile << PP.T_p[i][j];
			Tfile << " ";
		}
	}
	xfile.close();
	yfile.close();
	zfile.close();
	rfile.close();

	Tfile.close();

	for (int i = 0; i < grid; ++i)
	{
		for (int j = 0; j < PP.particle_count; ++j)
		{
			for (int k = 0; k < 15; ++k)
			{
				Neighfile << PP.neighbors[i][j][k];
				Neighfile << " ";
			}
			Neighfile << endl;
		}
		Neighfile << endl << endl;
	}
	Neighfile.close();
   // for (int u=0; u<PP.particle_count; u++)
   // {
   //   cout << PP.r_p[u] << " ";
   // }
}

#include<bits/stdc++.h>
#include "RT.h"
#include <fstream>

using namespace std;

// Function to read powder bed details from files
void PowderBed::GenerateBed(int grid_count, int particle_count)
{
	ifstream xfile("../ThAn/Results/x.txt");
	ifstream yfile("../ThAn/Results/y.txt");
	ifstream zfile("../ThAn/Results/z.txt");
	ifstream rfile("../ThAn/Results/r.txt");

	ifstream Neighfile("../ThAn/Results/Neighbors.txt");

	for (int i = 0; i < PP.particle_count; ++i)
	{
		rfile >> PP.r_p[i];
	}
	for (int i = 1; i <= PP.grid_count; ++i)
	{
		for (int j = 0; j < PP.particle_count; ++j)
		{
			xfile >> PP.x_p[i][j];
			yfile >> PP.y_p[i][j];
			zfile >> PP.z_p[i][j];
		}
	}
	xfile.close();
	yfile.close();
	zfile.close();
	rfile.close();

	for (int i = 1; i <= PP.grid_count; ++i)
	{
		for (int j = 0; j < PP.particle_count; ++j)
		{
			for (int k = 0; k < 15; ++k)
			{
				Neighfile >> PP.neighbors[i][j][k];
			}
		}
	}
	Neighfile.close();
}
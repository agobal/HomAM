#include<bits/stdc++.h>
#include "RT.h"
#include <fstream>

using namespace std;

// Function to read powder bed details from files
void PowderBed::GenerateBed(int grid_count, int particle_count)
{
	ifstream xfile("./x.txt");
	ifstream yfile("./y.txt");
	ifstream zfile("./z.txt");
	ifstream rfile("./r.txt");

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
}
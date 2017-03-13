#include<bits/stdc++.h>
#include "RT.h"
#include <ctime>
#include <fstream>
#include <iostream>

using namespace std;
/*
Ray tracing using DEM
// */

int main() 
{
	std::clock_t start;
	start = std::clock();
	// First, the program needs to read the x, y, z locations of powder particles, their respective radii r,
	// list of neighbors of particles and store them all in class arrays.
	// To create the desired array, we first need to determine array sizes. A function does that by reading
	// the array size from the x file.

	// PowderBed pb;
	int grid_count = 0;
	int particle_count = 0;
    string line;
 
    /* Creating input filestream */ 
    ifstream file("../ThAn/Results/x.txt");
    while (getline(file, line))
    {
    	if (grid_count == 0)
    	{
			particle_count = std::distance(
							 	std::istream_iterator<std::string>(std::istringstream(line) >> std::ws),
		                    	std::istream_iterator<std::string>());
    	}
        grid_count++;
    }

    cout << "Number of grids:" << grid_count << "\t" << "Number of particles:" << particle_count << "\n";

    cout << "Size of powder bed determined!" << "\t" << "time: " << ( std::clock() - start ) / (double)(CLOCKS_PER_SEC/1000) << "ms" << endl;

	// Constructor for the PowderBed class
	PowderBed pb(grid_count, particle_count);

	cout << "Required arrays constructed!" << "\t" << "time: " << ( std::clock() - start ) / (double)(CLOCKS_PER_SEC/1000) << "ms" << endl;

	// Next we need to read the x, y, z, r and neighbor values into the dynamic arrays
	pb.GenerateBed(grid_count, particle_count);

	// Perform the ray tracing analysis
	pb.RayTracing();

    return 0;
}
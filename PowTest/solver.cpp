#include<bits/stdc++.h>
#include "test.h"
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
	int grid_count = 100;
	int particle_count = 1;
    string line;

    cout << "Number of grids:" << grid_count << "\t" << "Number of particles:" << particle_count << "\n";

    cout << "Size of powder bed determined!" << "\t" << "time: " << ( std::clock() - start ) / (double)(CLOCKS_PER_SEC/1000) << "ms" << endl;

	// Constructor for the PowderBed class
	PowderBed pb(grid_count, particle_count);

	cout << "Required arrays constructed!" << "\t" << "time: " << ( std::clock() - start ) / (double)(CLOCKS_PER_SEC/1000) << "ms" << endl;

	// Next we need to read the x, y, z, r and neighbor values into the dynamic arrays
	pb.GenerateBed(grid_count, particle_count);
	cout << grid_count << "\t" << particle_count << "\n";
	cout << "Bed Generated!" << "\t" << "time: " << ( std::clock() - start ) / (double)(CLOCKS_PER_SEC/1000) << "ms" << endl;
	pb.Dele();
    return 0;
}
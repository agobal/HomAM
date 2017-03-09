#include<bits/stdc++.h>
#include "TrussDEM.h"
#include <ctime>
#include <fstream>
#include <iostream>

using namespace std;
/*
This program solves the 3d elasticity problem on the powder bed using a discrete element
method that considers each powder particle as a node in a 3d truss. The stiffness and mass matrices of the 
truss elements are calculated using results of current literature
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

    cout<<"time: " << ( std::clock() - start ) / (double)(CLOCKS_PER_SEC/1000) << "ms" << endl;

	// Constructor for the PowderBed class
	PowderBed pb(grid_count, particle_count);

	cout<<"time: " << ( std::clock() - start ) / (double)(CLOCKS_PER_SEC/1000) << "ms" << endl;
    return 0;
}
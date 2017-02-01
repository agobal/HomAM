#include <fstream>
#include <iostream>
#include <stdio.h>
#include "SLS.h"
#include <ctime>
#include <math.h>
#include <typeinfo>
#include <vector>
#include <algorithm>
#include <cstring>

using namespace std;

/* 
This program tries to perform a 3D full scale simulation of the selective laser sintering process.
A multiscale adaptive discrete element method is used for performing thermal simulations and the laser 
beam is modeled as a Gaussian heat input on top of the powder bed.
*/

int main()
{	
	/* Starting the chronometer */
	int start_s = clock();

	BedGeometry BG;
	/* Powder bed dimensions and layer thickness (all in milimeters) */

	BG.bed_x = 1.0;
	BG.bed_y = 1.0;
	BG.bed_z = 2.0;
	BG.layer_thickness = 0.05;
	BG.layer_thickness_big = 0.0625;

	/* Grid size and number of grids (all in milimeters) */

	BG.grid_x = 0.1;	// Size of grid-x
	BG.grid_y = 0.1;	// Size of grid-y
	BG.grid_z = 0.05;	// Size of grid-z
	BG.grid_z_big = 0.0625;

	BG.grid_volume = BG.grid_x*BG.grid_y*BG.grid_z;		// Volume of each grid
	BG.grid_volume_big = BG.grid_x*BG.grid_y*BG.grid_z_big;		// Volume of each grid

	BG.num_grid_x = int(BG.bed_x/BG.grid_x);		// Number of grids in x-direction
	BG.num_grid_y = int(BG.bed_y/BG.grid_y);		// Number of grids in y-direction
	BG.num_grid_z = int(BG.layer_thickness/BG.grid_z);		// Number of grids in z-direction
	BG.num_grid_z_big = int(BG.layer_thickness_big/BG.grid_z_big);		// Number of grids in z-direction

	static int num_grid = BG.num_grid_x*BG.num_grid_y*BG.num_grid_z;	// Total number of grids in the bed
	static int num_grid_big = BG.num_grid_x*BG.num_grid_y*BG.num_grid_z_big;	// Total number of grids in the bed
	BG.num_grid = num_grid;
	BG.num_grid_big = num_grid_big;

	BG.bed_volume = BG.bed_x*BG.bed_y*BG.bed_z;

	//////////////////////////////////////////////////////
	//////////// Packing of the particle level ///////////
	//////////////////////////////////////////////////////

	/* Powder packing properties */
	ParticleChar PC1;	// Loading the powder packing data structre from SLS.h

	PC1.avgrd = 7.5/1000.0;		// Average diameter of particles inside the packing
	PC1.stddev = 2.5/1000.0;		// Standard deviation of the particles inside packing
	PC1.packfrac = 0.63;	// Packing fraction of the bed
	PC1.nmin = BG.grid_volume/((4.0/3.0)*4.0*atan(1.0)*pow((PC1.avgrd), 3.0));		// Minimum number of particles inside each grid

	// Finding out the number of powders inside each grid 
	float* p_r;
	p_r = RadiusFinder(PC1, BG.grid_volume);

	static int num_particle_grid = int(p_r[0]);

	// Creating the new array for particle diameters
	float particle_radius[num_particle_grid];
	for (int i = 0; i < num_particle_grid; ++i)
	{
		particle_radius[i] = p_r[i + 1];
	}
	// Set the powder bed structure lol
	PowderBed PB;
	PB.particle_count = num_particle_grid;
	PB.cell_count = BG.num_grid;

	// Copy the particle radii to the class
	memcpy(PB.r_particles, particle_radius, sizeof(particle_radius));

	// Populate the powder bed properties class
	PB = PackingGenerator(PC1, BG, PB);
	cout << PB.cell_count << " " << PB.particle_count << endl;

	//////////////////////////////////////////////////////
	//////////// Packing of the element level ////////////
	//////////////////////////////////////////////////////

	/* Powder packing properties */
	ParticleChar PC_BIG;	// Loading the powder packing data structre from SLS.h

	PC_BIG.avgrd = 20/1000.0;		// Average diameter of particles inside the packing
	PC_BIG.stddev = 5.0/1000.0;		// Standard deviation of the particles inside packing
	PC_BIG.packfrac = 0.63;	// Packing fraction of the bed
	PC_BIG.nmin = BG.grid_volume_big/((4.0/3.0)*4.0*atan(1.0)*pow((PC_BIG.avgrd), 3.0));		// Minimum number of particles inside each grid

	/* Finding out the number of powders inside each grid */
	float* p_r_b;
	p_r_b = RadiusFinder(PC_BIG, BG.grid_volume_big);

	static int num_big_particle_grid = int(p_r_b[0]);

	// Creating the new array for particle diameters
	float big_particle_radius[num_big_particle_grid];
	for (int i = 0; i < num_big_particle_grid; ++i)
	{
		big_particle_radius[i] = p_r_b[i + 1];
	}
	// Set the powder bed structure lol
	PowderBed PB_BIG;
	PB_BIG.particle_count = num_big_particle_grid;
	PB_BIG.cell_count = BG.num_grid_big;

	// Copy the particle radii to the class
	memcpy(PB_BIG.r_particles, big_particle_radius, sizeof(big_particle_radius));

	// Populate the powder bed properties class
	PB_BIG = PackingGenerator(PC_BIG, BG, PB_BIG);

	//////////////////////////////////////////////////////
	/////////////// Heat transfer analysis ///////////////
	//////////////////////////////////////////////////////

	// Function for heat transfer analysis
	int output_timestep = 500;
	TempProfile TP;
	TP = LaserSintering(PB, PB_BIG, output_timestep);

	// Write temperature values to text file
	ofstream myfile;
	myfile.open("temp.txt");
	for (int cell = 0; cell < PB.cell_count; ++cell)
	{
		for (int p = 0; p < PB.particle_count; ++p)
		{
			myfile << TP.T[cell][p] << " ";
		}
		myfile << endl;
	}
	myfile.close();
	/* Making up the powder bed */
	// PB = PackingGenerator(PC1, particle_diameter, );	// Generating the desired packing 

/*	//Making up the smaller packing within the large one
	SPowderBed SBed;
	SBed = SmallPackingGenerator(Bed);

	cout<<"Number of particles:"<<SBed.cnt<<endl;

	//Generating the shape function matrices and element stiffness and capacity matrices
	ElementCoefficients EC;
	EC = ShFcnGen(Bed, SBed);

	//Finding elements on sides



	//Writing powder bed data into .txt files
	OutputWriter(Bed);

	int stop_s = clock();
	cout<<"Stiffness matrices created"<<endl;*/
	// cout << "time: " << (stop_s - start_s)/double(CLOCKS_PER_SEC) << endl;

	return 0;
}
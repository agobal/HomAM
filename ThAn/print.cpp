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
	ofstream xfile("./Results/x.txt");
	ofstream yfile("./Results/y.txt");
	ofstream zfile("./Results/z.txt");
	ofstream rfile("./Results/r.txt");

	ofstream Tfile("./Results/T.txt");

	ofstream Neighfile("./Results/Neighbors.txt");

	for (int i = 0; i < PP.particle_count; ++i)
	{
		rfile << PP.r_p[i];
		rfile << " ";
	}
	for (int i = 1; i <= PP.grid_count; ++i)
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
		xfile << "\n";
		yfile << "\n";
		zfile << "\n";
		Tfile << "\n";
	}
	xfile.close();
	yfile.close();
	zfile.close();
	rfile.close();

	Tfile.close();

	for (int i = 1; i <= grid; ++i)
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

// Save location of the large particle
	ofstream Xfile("./Results/Xl.txt");
	ofstream Yfile("./Results/Yl.txt");
	ofstream Zfile("./Results/Zl.txt");
	ofstream Rfile("./Results/Rl.txt");

	ofstream Tlfile("./Results/Tl.txt");

	ofstream Neighfilel("./Results/Neighborsl.txt");

	for (int i = 0; i < PPl.particle_count; ++i)
	{
		Rfile << PPl.r_p[i];
		Rfile << " ";
	}
	for (int i = 1; i <= PPl.grid_count; ++i)
	{
		for (int j = 0; j < PPl.particle_count; ++j)
		{
			Xfile << PPl.x_p[i][j];
			Xfile << " ";
			Yfile << PPl.y_p[i][j];
			Yfile << " ";
			Zfile << PPl.z_p[i][j];
			Zfile << " ";
			Tlfile << PPl.T_p[i][j];
			Tlfile << " ";
		}
		Xfile << "\n";
		Yfile << "\n";
		Zfile << "\n";
		Tlfile << "\n";
	}
	Xfile.close();
	Yfile.close();
	Zfile.close();
	Rfile.close();

	Tfile.close();

	for (int i = 1; i <= grid; ++i)
	{
		for (int j = 0; j < PPl.particle_count; ++j)
		{
			for (int k = 0; k < 15; ++k)
			{
				Neighfilel << PPl.neighbors[i][j][k];
				Neighfilel << " ";
			}
			Neighfilel << endl;
		}
		Neighfilel << endl << endl;
	}
	Neighfilel.close();
}

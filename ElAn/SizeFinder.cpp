#include<bits/stdc++.h>
#include <iostream>
#include <fstream>
#include "TrussDEM.h"

using namespace std;

void SizeFinder(int grid_count, int particle_count)
{
    // int grid_count = 0;
    // int particle_count = 0;
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
}
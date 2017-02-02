#include<bits/stdc++.h>
#include "LS.h"

using namespace std;
 

 
PowBed::PowBed(int particle_count, int grid_count)
{
    this->particle_count = particle_count;
    this->grid_count = grid_count;
 
    // Create a dynamic array of pointers
    x_p = new float* [grid_count];
    y_p = new float* [grid_count];
    z_p = new float* [grid_count];

    // Create a row for every pointer
    for (int i=0; i<grid_count; i++)
    {
       // Note : Rows may not be contiguous
       x_p[i] = new float[particle_count];
       y_p[i] = new float[particle_count];
       z_p[i] = new float[particle_count];
 
       // Initialize all locations are zero
       memset(x_p[i], 0.0, particle_count*sizeof(float));
       memset(y_p[i], 0.0, particle_count*sizeof(float));
       memset(z_p[i], 0.0, particle_count*sizeof(float));
    }

    r_p = new float [particle_count];
    memset(r_p, 0.0, particle_count*sizeof(float));
}
 
// Utility method to print adjacency matrix
void PowBed::print()
{
   for (int u=0; u<grid_count; u++)
   {
      for (int v=0; v<particle_count; v++)
         cout << x_p[u][v] << " ";
      cout << endl;
   }
}

// Set initial locations
void PowBed::SetInitials(int grid_no, int particle_no)
{
	x_p[grid_no][particle_no] = 15;
}
 
// Driver method
int main()
{
    // Create a graph given in the above diagram
    PowBed g(500, 10000);
    g.SetInitials(2,4);
 
    g.print();
 
    return 0;
}
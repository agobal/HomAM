#include<bits/stdc++.h>
#include "LS.h"

using namespace std;
 
// Driver method
int main()
{
    // Create a graph given in the above diagram
    // PowBed::ParProps PP;
    int par = 100;
    int grid = 128;
    int par_l = 15;

    PowBed g(grid, par, par_l);


    g.SetInitials(grid, par, par_l);

    float power = 2.0; // Laser power (watts)
    float speed = 0.001; // Laser beam speed (m/s)
    g.Sintering(power, speed);

    g.print();
 	g.PowDel();
    return 0;
}
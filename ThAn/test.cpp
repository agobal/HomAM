#include<bits/stdc++.h>
#include "LS.h"

using namespace std;
 
// Driver method
int main()
{
    // Create a graph given in the above diagram
    // PowBed::ParProps PP;
    int par = 100;
    int grid = 64;
    PowBed g(grid, par);


    g.SetInitials(grid, par);

    float power = 50; // Laser power (watts)
    float speed = 1.5; // Laser beam speed (m/s)
    g.Sintering(power, speed);

    g.print();
 	g.PowDel();
    return 0;
}
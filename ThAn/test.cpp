#include<bits/stdc++.h>
#include "LS.h"
#include <ctime>

using namespace std;
 
// Driver method
int main()
{
    // Create a graph given in the above diagram
    // PowBed::ParProps PP;
    int par = 100;
    int grid = 64;
    int par_l = 15;

    PowBed g(grid, par, par_l);

  std::clock_t start;
  start = std::clock();
    g.SetInitials(grid, par, par_l);
  cout<<"time: " << ( std::clock() - start ) / (double)(CLOCKS_PER_SEC/1000) << "ms" << endl;




    float power = 25; // Laser power (watts)
    float speed = 1; // Laser beam speed (m/s)

  std::clock_t start;
  start = std::clock();
    g.Sintering(power, speed);
  cout<<"time: " << ( std::clock() - start ) / (double)(CLOCKS_PER_SEC/1000) << "ms" << endl;

    g.print();
 	g.PowDel();
    return 0;
}
#include<bits/stdc++.h>
#include "LS.h"

using namespace std;
 
// Driver method
int main()
{
    // Create a graph given in the above diagram
    // PowBed::ParProps PP;
    int par = 100;
    int grid = 3375;
    PowBed g(grid, par);


    g.SetInitials(grid, par);
    g.print();
 	g.PowDel();
    return 0;
}
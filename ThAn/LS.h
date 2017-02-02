// A Class to represent the powder bed arrangement
class PowBed
{
    int particle_count;    // No. of particles
 	int grid_count;			// No. of grids

    // arrays for x, y and z location of the powder particles
    float **x_p;
    float **y_p;
    float **z_p;
    float *r_p;
 
public:
    PowBed(int particle_count, int grid_count);   // Constructor
 
    // function to add an edge to graph
    void SetInitials(int grid_no, int particle_no);
    void print();
};
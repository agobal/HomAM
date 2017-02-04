// A Class to represent the powder bed arrangement
class PowBed
{
public:
    int grid;
    int par;
    int num_grid_x;
    int num_grid_y;
    int num_grid_z;
    float grid_x;
    float grid_y;
    float grid_z;
    // arrays for x, y and z location of the powder particles
    struct ParProps
    {
        int particle_count;    // No. of particles
        int grid_count;         // No. of grids
        float avgrd;            // Average particle radius
        float stddev;           // Average standard deviation of particle radii
        float packfrac;         // Packing fraction of the powder bed
        float **x_p;
        float **y_p;
        float **z_p;
        int **sintering_flag;
        float *r_p;
    };
    ParProps PP;

    PowBed(int grid, int par);   // Constructor

    // function to add an edge to graph
    void SetInitials(int grid, int par);
    void print();
    void PowDel();   // Delete
};
// Function to determine array sizes
// void SizeFinder(int grid_count, int particle_count);

class PowderBed
{
// Bulk properties
public:
	int grid_count;
	int particle_count;
    // Bulk properties of the powder bed
    int grid;
    int par;
    int par_l;
    int num_grid_x;
    int num_grid_y;
    int num_grid_z;
    float grid_x;
    float grid_y;
    float grid_z;
    // arrays for x, y and z location of the powder particles
    struct PowProps
    {
        float avgrd;            // Average particle radius
        float stddev;           // Average standard deviation of particle radii
        float packfrac;         // Packing fraction of the powder bed
        int particle_count;    // No. of particles
        int grid_count;         // No. of grids
        float **x_p;
        float **y_p;
        float **z_p;
        float *r_p;
    };
    PowProps PP;

    PowderBed(int grid_count, int particle_count);   // Constructor

    void GenerateBed(int grid_count, int particle_count);

    void Dele();
};
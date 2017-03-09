// Function to determine array sizes
// void SizeFinder(int grid_count, int particle_count);

class PowderBed
{
// Bulk properties
public:
	int grid_count;
	int particle_count;

    // arrays for x, y and z location of the powder particles
    struct PowProps
    {
        int particle_count;    // No. of particles
        int grid_count;         // No. of grids
        float **x_p;
        float **y_p;
        float **z_p;
        int ***neighbors;
        float *r_p;

        // Properties for elasticity simulation
        float *u;	// Displacement dof (size particle_c*grid_c*dof)
        float *f;  // Force acted upon particles (size particle_c*grid_c*dof)
        float **K; // Stiffness matrix of the system
        float **M; // Mass matrix of the system
    };
    PowProps PP;

    PowderBed(int grid_count, int particle_count);   // Constructor
};
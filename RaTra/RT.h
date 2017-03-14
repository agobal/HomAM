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
        float *r_p;
    };
    PowProps PP;

    struct LaserRays
    {
        int num_rays;
        float *ray_power;
        float **intersection_point;
    };
    LaserRays LR;

    PowderBed(int grid_count, int particle_count);   // Constructor

    void GenerateBed(int grid_count, int particle_count);

    void RayTracing();

    void Dele();
};
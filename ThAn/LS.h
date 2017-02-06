// A Class to represent the powder bed arrangement
class PowBed
{
public:
    // Bulk properties of the powder bed
    int grid;
    int par;
    int num_grid_x;
    int num_grid_y;
    int num_grid_z;
    float grid_x;
    float grid_y;
    float grid_z;

    // Mechanical properties of powders
    float E; // Young's modulus
    float nu; // Poisson's ratio
    float k1; // Thermal conductivity of particle 1
    float k2; // Thermal conductivity of particle 2
    float k_air; // Thermal conductivity of air

    // Sintering properties
    float W;    // Laser power
    float dt;   // Sintering timestep

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
        int ***neighbors;
        float *r_p;
        float **T_p;            // Particle temperature
        float **delta_T_p;      // delta particle temperature (over time)
    };
    ParProps PP;

    PowBed(int grid, int par);   // Constructor

    // function to set initial location of powder particles
    void SetInitials(int grid, int par);

    // Function for carrying out the heat transfer simulation
    void Sintering(float powder, float speed);

    // Function to print results onto the screen
    void print();
    // Function to clear the heap memory after the program's done
    void PowDel();   // Delete
};
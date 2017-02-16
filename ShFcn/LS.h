// A Class to represent the powder bed arrangement
class PowBed
{
public:
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

    // Sintering properties
    float W;    // Laser power
    float dt;   // Sintering timestep
    float K_F;
    float Laser_Intensity; 

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
        float **T_temp;      // delta particle temperature (over time)
        float **E;          // Internal energy of powder particles
        int *cell_sinter_flag; // For fast computations
        float *max_temp;
        float **K_F;
    };
    ParProps PP;
    ParProps PPl;

    PowBed(int grid, int par, int par_l);   // Constructor

    // function to set initial location of powder particles
    void SetInitials(int grid, int par, int par_l);

    // Function for carrying out the heat transfer simulation
    void Sintering(float powder, float speed);

    // Function for calculating the heat transfer coefficient between powders and their location change
    void CondCoeff(int cell1, int par1, int cell2, int par2, float delta_t);

    // Function to print results onto the screen
    void print();
    // Function to clear the heap memory after the program's done
    void PowDel();   // Delete
};
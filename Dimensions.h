#ifndef DIMENSIONS
#define DIMENSIONS


// These two parameters define the h value in the CFL condition of del t = h /u(max)
// h in example below is 0.2 [coz 8m/40 = 20cm = 0.2m (velocity is in m/s)]
#define PHYSICAL_WIDTH 1
#define PHYSICAL_HEIGHT 1

#define GRID_NI 40
#define GRID_NJ 40

// this defines the cell size (used for drawing purposes chiefly)
#define RENDER_SIZE_CELL 20

// the fuel and hot gas densities [unit is kg/(m*m*m)]
#define DENSITY_AIR   0.0
#define DENSITY_FUEL  1.0
#define DENSITY_HOTGAS 0.01


#define FLAME_HEIGHT 0.5

//for testing purposes only
#define V_MAX 70
#define U_MAX 70

// the alpha parameter [unit is meters/Kelvin*second*second]
#define TEMP_ALPHA 0.0075
//[unit of CT is K/s]
#define CT 3000 


//auxillary
#define FULL_SCREEN	false

#define START_TIME      0.0
#define END_TIME        10.0



#endif


#ifndef SIM_H_
#define SIM_H_

extern "C" {
	#include "config.h"
}


// These two parameters define the h value in the CFL condition of del t = h /u(max)
// h in example below is 0.2 [coz 8m/40 = 20cm = 0.2m (velocity is in m/s)]
#define PHYSICAL_WIDTH 		config_get_int("PHYSICAL_WIDTH")
#define PHYSICAL_HEIGHT 	config_get_int("PHYSICAL_HEIGHT")

#define GRID_NI 	config_get_int("GRID_NI")
#define GRID_NJ 	config_get_int("GRID_NJ")

// this defines the cell size (used for drawing purposes chiefly)
#define RENDER_SIZE_CELL config_get_int("RENDER_SIZE_CELL")

// the fuel and hot gas densities [unit is kg/(m*m*m)]
#define DENSITY_AIR   	config_get_double("DENSITY_AIR)")
#define DENSITY_FUEL  	config_get_double("DENSITY_FUEL")
#define DENSITY_HOTGAS config_get_double("DENSITY_HOTGAS")


#define FLAME_HEIGHT config_get_double("FLAME_HEIGHT")

//for testing purposes only
#define V_MAX config_get_int("V_MAX")
#define U_MAX config_get_int("U_MAX")

// the alpha parameter [unit is meters/Kelvin*second*second]
#define TEMP_ALPHA 0.0075
//[unit of CT is K/s]
#define CT 3000


//auxillary
#define FULL_SCREEN	 false

#define START_TIME      0.0
#define END_TIME        10.0



#endif /* SIM_H_ */

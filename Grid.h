// $Id$
/**
 * @file Grid.hpp
 * Implements the staggered grid needed for the simulation.
 *
 * @brief Fluid Simulation -- Grid Declarations
 * @author Gaurav Bhagwat
 * @version 1.0
 */
// $Log$


#ifndef GRID_H_
#define GRID_H_

#include <boost/numeric/ublas/matrix.hpp>
#include "Vec.h"
#include "MetaConfig.h"
#include "Particles.h"

#include <vector>
#include <list>
#include <iostream>
using namespace std;

/**
 * Class Grid
 *
 * Implements the Grid for a fluid. Measurements like velocity, pressure, "phi" and/or Marker Particles
 * are the core components.
 *
 * Divergence, extrapolated velocity field, extrapolated phi (optional) and sparse matrices are auxillary data fields
 * needed.
 *
 * Initialization is a compulsory function that must be implemented to assign valid grid values before simulation.
 * The constructor must be the function to do so.
 */
using namespace boost::numeric::ublas;

class Grid {
public:

    // DATA MEMBERS

	bool makeFlammableBoundary;

    /// @brief the physical dimensions of the simulation domain (x-axis)
    double phyWidth;
    /// @brief the physical dimensions of the simulation domain (x-axis)
    double phyHeight;
    /// @brief number of cells (x-axis) for discretization j as columns horizontally
    int nj;
    /// @brief number of cells (y-axis) for discretization i as rows that go up vertically
    int ni;
    /// @brief cell_width =  phyWidth/ni
    double dx;

    /// @brief velocity component sampled in x-direction on staggered grid, size (ni+1)*nj
    matrix<double> u;

    /// @brief velocity component sampled in y-direction on staggered grid, size ni*(nj+1)
    matrix<double> v;
    
    /// @brief temp velocities may be helpful in advection
    matrix<double> temp_u, temp_v;
    
    matrix<Vec2d> velCenter;
    matrix<Vec2d> normal;
    matrix<double> omega;

    /// @brief advection weights
    matrix<double> advectionWeights;

    /// @brief marker cells
    matrix<int> marker;
    matrix<bool> flammable;
	
	/// @brief boundary cells
	matrix<int> boundary;

	///@brief levelSetPhi
	matrix<double> levelSetPhi;
	matrix<double> levelSetTemp;
	matrix<double> smokeDensity;
    

    /// @brief temperature for the flame
    matrix<double> temperature;
    

    /// @brief Vec2d(A double array of size=2) storing maxVelocity components maxVelocity[0] => u, maxVelocity[1] => v
    Vec2d maxVelocity;

    /// @brief pressure calculated during projection method, sampled at center of the grid cells
    matrix<double> pressure;
    matrix<double> temp_pressure;

    /// @brief (better known as divergence*some_constant, refer to Bridson's book) used in projection.
    matrix<double> rhs, customRhs, rhs_u_comp, rhs_v_comp;
    matrix<double> diagnose_u_rhs, diagnose_v_rhs;
    
    

    /// @brief auxillary matrices needed for pressure solve
    matrix<double>  Adiag, Aplusi, Aplusj;
    //matrix<double> myAdiag, Aplusy, Aplusx ;
    matrix<double> s, z;
    matrix <double> preConditioner, m;
    
    //matrix<double> myPreConditioner, q;
    //matrix<double> residualR, auxillaryZ, searchS;
    
    
    /// @brief particles for FLIP method
    //std::vector<Vec2d> flipParticles;

    /// @brief particles for capturing the fluid entities
     
    std::vector<Particles*> fireParticles;

    /// @brief radius of marker particles
    double particleRadius;

    /**
     * Constructor
     * @param ni
     * @param nj
     * @param phyWidth
     * @param phyHeight
     */
    Grid(int ni, int nj, double phyWidth, double phyHeight) {
        if (META_LOG)
            std::cout << "--META--\tConstructor\tGrid" << std::endl;

        this->dx = 0.0;
        this->particleRadius = 0.0;
        this->ni = ni;
        this->nj = nj;
        this->phyWidth = phyWidth;
        this->phyHeight = phyHeight;

        if (DEBUG_LOG)
            std::cout << "--DEBUG--\tGridConstructor\tni = " << ni << " nj = " << nj << std::endl;

    }

    /**
     * Destructor
     */
    virtual ~Grid() {
        if (META_LOG)
            std::cout << "--META--\tDestructor\tGrid" << std::endl;
    }

    /**
     * from base Class (need to be overridden)
     */
    void initializeGrid();

    /**
     * getter function
     * @param position
     * @return 
     */
    Vec2d getVelocity(const Vec2d& position);
    
    double getLevelSetPhi(const Vec2d& position);

    Vec2d trace_rk2(const Vec2d& position, double dt);
        
    


};

#endif /* GRID_H_ */
        

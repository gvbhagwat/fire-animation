#ifndef FLUIDSIM_H_
#define FLUIDSIM_H_

#include "Vec.h"
#include "Grid.h"

#include "MetaConfig.h"

#include <vector>
#include <iostream>

enum FluidBodyInitialization {
    STATIC_BED,
    DAM_BREAK,
    DOUBLE_DAM_BREAK,
    DAM_BREAK_OBSTACLE,
    CHANNEL_FLOW,
    CHANNEL_FLOW_OBSTACLE
};

class FluidSim {
protected:

    // Data Member
    Grid& rGrid;

public:

    FluidSim(Grid& grid) : rGrid(grid) {
        if (META_LOG)
            std::cout << "--META--\tConstructor\tFluidSim" << std::endl;
        //initialize(startTime, endTime, rGrid);
    }

    ~FluidSim() {
        if (META_LOG)
            std::cout << "--META--\tDestructor\tFluidSim" << std::endl;

    }


    /// @brief intialization routine, initialize time, fluid enity, solid boundary
    //    void initialize(Grid* grid, double startTime, double endTime);
    void initialize(double startTime, double endTime);

    /// @brief the core function "advancing" the simulation
    void advance(double timestep);

    /// @brief define the solid boundary
    void initializeSolidBoundaries();

    /// @brief initialize the fluid body as per required configurations, currently done by particles
    void initializeFluidBody(int choice);

    /// @brief the startTime and endTime of the simulation
    double startTime, endTime;

    /// @brief the time sub step fpr which the simulation is performed, enforced by cfl()
    double dt;

protected:

    /**
     * Advection routine.. encapsulates algorithm FLIP/PIC or simple
     * semi-Lagrangian advection
     * @param grid
     */
    void advect();

    /**
     * Add gravitational or buoyancy forces through out the domain
     * depending upon type of fluid
     * @param grid
     */
    void addForces();

    /**
     * Make fluids incompressible.
     * Modified Incomplete Cholesky Preconditioned Conjugate Gradient Algorithm.
     * @param grid
     */
    void project();

    void myProject();


    /**
     * The CFL condition calculates for minimum timestep
     * @param grid
     * @return timestep
     */
    double cfl();

    /**
     * PARTICLE BASED FLUID TRACING
     * @param grid
     */
    void advectParticles();

    /**
     * The marker and cell method's marker routines
     * @param grid
     */
    void markFluidCells();


    void applyBoundaryConditions();

    /**
     * @param grid
     */
    void markSolidCells();
    void clearMarkedFluidCells();



    // VELOCITY EXTRAPOLATION
    void constrainVelocity();


    // HELPER FUNCTIONS

    void calculateNegativeDivergence();
    void formCoefficientMatrixA();
    void formMyPreConditioner();
    void applyMyPreconditioner();
    void applyMyA();
    void solveForPressure();
    void updateVelocitiesFromPressure();
    
    // helper functions
    double dotProduct(matrix<double> m1, matrix<double> m2);

    void findDivergence();
    void formA();
    void applyA();
    void copyztos();
    double zdots();
    double zdotr();
    void formPreconditioner();
    void applyPreconditioner();
    double findRInform();
    void formRInform();
    void solvePressure();
    void applyPressure();


    void computeWeights();





};

#endif /* FLUIDSIM_H_ */


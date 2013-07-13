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

        //invalid value at initialization.. to be correctly given by the initialize()
        this->startTime = this->endTime = this->dt = -1.0;fbuoy = 9.8;

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

    void constructLevelSetPhi();

    /// @brief the startTime and endTime of the simulation
    double startTime, endTime;

    /// @brief the time sub step fpr which the simulation is performed, enforced by cfl()
    double dt;

    double fbuoy;

protected:

    /**
     * Advection routine.. encapsulates algorithm FLIP/PIC or simple
     * semi-Lagrangian advection
     *
     */
    void advect();

    /**
     * Add gravitational or buoyancy forces through out the domain
     * depending upon type of fluid
     *
     */
    void addForces();

    /**
     * @brief Make fluids incompressible.
     * Modified Incomplete Cholesky Preconditioned Conjugate Gradient Algorithm.
     *
     */
    void project();



    /**
     * @brief The CFL condition calculates for minimum timestep
     * @return timestep
     */
    double cfl();

    /**
     * PARTICLE BASED FLUID TRACING
     *
     */
    void advectParticles();

    void advectLevelSetPhi();

    /**
     * @brief The marker and cell method's marker routines
     *
     */
    void markFluidCells();

    /**
     * @brief construct level sets
     *
     */


    void applyBoundaryConditions();

    /**
     *
     */
    void markSolidCells();

    /**
     *
     */
    void clearMarkedFluidCells();





    // VELOCITY EXTRAPOLATION
    void constrainVelocity();

    void calculateNormal();
    void calculateCenterVelocityField();
    void confineVorticity();





//    void myProject();
//    // HELPER FUNCTIONS
//    void calculateNegativeDivergence();
//    void formCoefficientMatrixA();
//    void formMyPreConditioner();
//    void applyMyPreconditioner();
//    void applyMyA();
//    void solveForPressure();
//    void updateVelocitiesFromPressure();
//
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

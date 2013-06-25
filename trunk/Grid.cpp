// $Id$
/**
 * @file Grid.cpp
 * Implements the staggered grid needed for the simulation
 *
 * @brief Fluid Simulation Function Definitions
 * @author Gaurav Bhagwat
 * @version 1.0
 */
// $Log$
//#include "Dimensions.h"
#include "Grid.h"
#include "Vec.h"
#include "Utils.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

/**
 * Initializes the grid - All the matrices are zero-initialized
 * The currVelocityUcomp has (N+1) x N matrix size, while currVelocityVcomp has
 * N x (N+1) size. Also, the grid is now NxN and thus temperature, density and
 * pressure matrices are of the same size
 *
 * @param ni
 * @param nj
 *
 */
void Grid::initializeGrid() {
    if (META_LOG)
        cout << "--META--\tFunctionCall\tGrid::initializeGrid called" << endl;

    // must also be dy = phyHeight but I have implemented symmetric matrices.
    dx = phyWidth / (double) ni;

    //TODO my hack ..a very bad one (but quick).. memories were overlapping .. 
    int extra = 0;

    // nj or nx number of cells (x-axis) for discretization j as columns horizontally
    // ni or ny number of cells (y-axis) for discretization j as columns horizontally

    // need to have an extra column horizontally or in x direction
    // but the same number of rows in y direction or ni
    // i will always be rows in matrix (that is y ) and j will always be column in matrix
    // that is to the right horizontal
    u.resize(ni + extra, nj + 1 + extra, false);
    temp_u.resize(ni + extra, nj + 1 + extra, false);


    // v needs extra row or extra information in y axis
    // hence the dimensions are ni+1 rows and the same number of columns
    // or nj 
    v.resize(ni + 1 + extra, nj + extra, false);
    temp_v.resize(ni + 1 + extra, nj + extra, false);


    // always calculated at cell centers hence only ni rows in y direction
    // and nj columns in x direction
    pressure.resize(ni + extra, nj + extra, false);
    temp_pressure.resize(ni + extra, nj + extra, false);

    // another information stored at the center of the cell hence ni x nj
    marker.resize(ni + extra, nj + extra);
	this->boundary.resize(ni + extra, nj +extra);
	this->boundary.resize(ni + extra, nj + extra);


    // myown test matrices



    // next come the pressure solve matrices
    rhs.resize(ni + extra, nj + extra);
    rhs_u_comp.resize(ni + extra, nj + extra);
    rhs_v_comp.resize(ni + extra, nj + extra);

    s.resize(ni + extra, nj + extra);
    z.resize(ni + extra, nj + extra);
    preConditioner.resize(ni + extra, nj + extra);
    m.resize(ni + extra, nj + extra);


    Adiag.resize(ni + extra, nj + extra);
    Aplusi.resize(ni + extra, nj + extra);
    Aplusj.resize(ni + extra, nj + extra);

    // for diagnostic purposes

    this->customRhs.resize(ni + extra, nj + extra);
    this->diagnose_u_rhs.resize(ni + extra, nj + extra);
    this->diagnose_v_rhs.resize(ni + extra, nj + extra);

    this->myAdiag.resize(ni + extra, nj + extra);
    this->Aplusy.resize(ni + extra, nj + extra);
    this->Aplusx.resize(ni + extra, nj + extra);

    this->myPreConditioner.resize(ni + extra, nj + extra);
    this->q.resize(ni + extra, nj + extra);
    this->residualR.resize(ni + extra, nj + extra);
    this->searchS.resize(ni + extra, nj + extra);
    this->auxillaryZ.resize(ni + extra, nj + extra);

    this->advectionWeights.resize(ni + 3 + extra, nj + 3 + extra);

    // particle radius for particle tracking and maxVelocity 
    // simple initialization
    particleRadius = dx / sqrt(2.0);
    maxVelocity[0] = maxVelocity[1] = 0.0;

    if (SIM_LOG)
        cout << "--SIMLOG--\tGrid Initialized" << endl;
}

/**
 * Interpolate velocity from the MAC grid.
 * @param position
 * @return Velocity vector Vec2d
 */
Vec2d Grid::getVelocity(const Vec2d& position) {
    if (META_LOG)
        std::cout << "--META--\tFunctionCall\tRenderer::getVelocity"
            << std::endl;

    int i, j;
    double fx, fy;

    Vec2d posX = position / dx - Vec2d(0.0, 0.5);
    Vec2d posY = position / dx - Vec2d(0.5, 0.0);

    // gets the i,j value for a given position and fraction value
    // gets the fx, fy to used in interpolation.

    get_barycentric(posX[0], i, fx, 0, ni);
    get_barycentric(posX[1], j, fy, 0, nj);

    //    this->advectionWeights(i,j) += fx;
    //    this->advectionWeights(i+1,j) += fx;
    //    this->advectionWeights(i,j+1) += fx;
    //    this->advectionWeights(i+1,j+1) += fx;
    //    
    //    cout<<"Initial Position ="<<position[0]<<" and "<<position[1]<<endl;
    //    cout<<"posX value = "<<posX[0]<<" and "<<posX[1]<<endl;
    //    cout<<"I_VALUE for u "<< i <<"J_VALUE for u"<< j<< endl;


    double u_value = bilerp(
            u(i, j), u(i + 1, j),
            u(i, j + 1), u(i + 1, j + 1),
            fx, fy);

    //    cout<<"u(i,j) --" << u(i,j)<<endl;
    //    cout<<"u(i+1,j) --" << u(i+1,j)<<endl;
    //    cout<<"u(i,j+1) --" << u(i,j+1)<<endl;
    //    cout<<"u(i+1,j+1) --" << u(i+1,j+1)<<endl;
    //    cout<<" fx = "<<fx<<" and "<<" fy = "<<fy<<endl;



    // This function also handles the boundary cases
    get_barycentric(posY[0], i, fx, 0, ni);
    get_barycentric(posY[1], j, fy, 0, nj);


    //    cout<<"posY value = "<<posY[0]<<" and "<<posY[1]<<endl;
    //    cout<<"I_VALUE for v "<< i <<"J_VALUE for v"<< j<< endl;
    //     
    //    
    //    
    //     cout<<"v(i,j) --" <<v(i,j)<<endl;
    //     cout<<"v(i+1,j) --" <<v(i+1,j)<<endl;
    //     cout<<"v(i,j+1) --" <<v(i,j+1)<<endl;
    //     cout<<"v(i+1,j+1) --" <<v(i+1,j+1)<<endl;
    //     cout<<" fx = "<<fx<<" and "<<" fy = "<<fy<<endl;
    if (fx < 1) {
        this->advectionWeights(i, j) += fx;
        this->advectionWeights(i + 1, j) += fx;
        this->advectionWeights(i, j + 1) += fx;
        this->advectionWeights(i + 1, j + 1) += fx;
    }

    double v_value = bilerp(
            v(i, j), v(i + 1, j),
            v(i, j + 1), v(i + 1, j + 1),
            fx, fy);

    return Vec2d(u_value, v_value);
}

Vec2d Grid::trace_rk2(const Vec2d& position, double dt) {
    Vec2d velocity = getVelocity(position);
    velocity = getVelocity(position + 0.5f * dt * velocity);
    Vec2d newPosition = position + dt*velocity;
    return newPosition;
}
                        

/* 
 * File:   GridTest.cpp
 * Author: gaurav
 * 
 * Created on 22 October, 2012, 6:03 PM
 */


#include "Utils.h"
#include "Vec.h"
#include "MetaConfig.h"
#include "sim.h"
#include "Grid.h"

#include <iostream>
#include <climits>
#include <cmath>
#include <cstdio>
using namespace std;

#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;

/**
 * 
 * @return 
 */
bool checkInitialize(Grid* rGrid) {

    cout << "Checking Grid Initialization" << endl;

    bool dimTest = GRID_NI == rGrid->ni && GRID_NJ == rGrid->nj;

    if (dimTest)
        cout << "PASS\tdimension assignment (ni,nj) " << "(" << rGrid->ni << "," << rGrid->nj << ")" << endl;
    else
        cout << "FAIL\tdimension assignment (ni,nj) " << "(" << rGrid->ni << "," << rGrid->nj << ")" << endl;

    double tdx = (double) PHYSICAL_WIDTH / (double) GRID_NI;
    bool dxTest = tdx == rGrid->dx;
    if (dxTest)
        cout << "PASS\tdx value (tdx,rGrid->dx) " << "(" << tdx << "," << rGrid->dx << ")" << endl;
    else
        cout << "FAIL\tdx value (tdx,rGrid->dx) " << "(" << tdx << "," << rGrid->dx << ")" << endl;


    unsigned int sizeUx = GRID_NI + 1;
    unsigned int sizeUy = GRID_NJ;

    bool uTest = sizeUx == rGrid->u.size1() && sizeUy == rGrid->u.size2();
    if (uTest)
        cout << "PASS\tu velocity dimensions (size1,size2) == " << "(" << rGrid->u.size1() << "," << rGrid->u.size2() << ")" << endl;
    else
        cout << "FAIL\tu velocity dimensions (" << sizeUx << ", " << sizeUy << ") != " << "(" << rGrid->u.size1() << "," << rGrid->u.size2() << ")" << endl;


    unsigned int sizeVx = GRID_NI;
    unsigned int sizeVy = GRID_NJ + 1;

    bool vTest = sizeVx == rGrid->v.size1() && sizeVy == rGrid->v.size2();
    if (vTest)
        cout << "PASS\tv velocity dimensions (size1,size2) == " << "(" << rGrid->v.size1() << "," << rGrid->v.size2() << ")" << endl;
    else
        cout << "FAIL\tv velocity dimensions (" << sizeVx << ", " << sizeVy << ") != " << "(" << rGrid->v.size1() << "," << rGrid->v.size2() << ")" << endl;


    return (vTest && uTest && dimTest && dxTest);

}

/**
 * 
 * @return 
 */
bool checkGetVelocity(Grid* rGrid) {
    // assigning test velocities

    int x, y;

    x = 0;
    y = 0;
    Vec2d testPos((x + 0.5) * rGrid->dx, (y + 0.5) * rGrid->dx);

    rGrid->u(x, y) = 10.0;
    rGrid->u(x + 1, y) = 10.0;
    rGrid->u(x + 1, y + 1) = 10.0;
    rGrid->u(x, y + 1) = 10.0;

    rGrid->v(x, y) = 20.0;
    rGrid->v(x + 1, y) = 20.0;
    rGrid->v(x + 1, y + 1) = 20.0;
    rGrid->v(x, y + 1) = 20.0;

    Vec2d result = rGrid->getVelocity(testPos);

    if (result[0] == 10 && result[1] == 20)
        cout << "PASS\tu_velocity\ttrivial test" << endl;
    else {
        cout << "FAIL\tu_velcoity\ttrivial test" << endl;
        cout << "\t\tFail result " << result[0] << " " << result[1] << endl;
    }

    x = 39;
    y = 39;
    testPos[0] = ((x + 0.5) * rGrid->dx);
    testPos[1] = ((y + 0.5) * rGrid->dx);

    rGrid->u(x, y) = 10.0;
    rGrid->u(x + 1, y) = 10.0;
    rGrid->u(x + 1, y + 1) = 10.0;
    rGrid->u(x, y + 1) = 10.0;

    rGrid->v(x, y) = 20.0;
    rGrid->v(x + 1, y) = 20.0;
    rGrid->v(x + 1, y + 1) = 20.0;
    rGrid->v(x, y + 1) = 20.0;

    result = rGrid->getVelocity(testPos);

    if (result[0] == 10 && result[1] == 20) {
        cout << "PASS\tu_velocity\tBOREDER test" << endl;
        cout << result << endl;
    } else {
        cout << "FAIL\tu_velcoity\tBORDER test" << endl;
        cout << "\t\tFail result " << result[0] << " " << result[1] << endl;
        cout << result << endl;
    }

    return false;
}

/**
 * 
 * @return 
 */
bool checkRenderingCoordinates(Grid* rGrid) {
    for (int j = 0; j < rGrid->nj; ++j) for (int i = 0; i < rGrid->ni; ++i) {

            Vec2d pos((i + 0.5) * rGrid->dx, (j + 0.5) * rGrid->dx);
            Vec2d vel = rGrid->getVelocity(pos);

            // per position case study starts here

            std::cout << "-----------------------" << std::endl;
            std::cout << "pos  = " << pos << " vel = " << vel;


            //vel = vel * 0.0005;
            if (mag(vel)) {
                vel = vel / mag(vel);
                vel = vel * rGrid->dx *
                		0.5;

            }
            vel = pos + vel;

            std::cout << " arrow = " << vel / rGrid->dx << std::endl;
            std::cout << "Actual pos (base) = " << pos / rGrid->dx;
            std::cout << "; Actual arrow (tip) = " << vel / rGrid->dx << std::endl;
            std::cout << "rendering (base) = " << pos / (rGrid->dx / RENDER_SIZE_CELL);
            std::cout << "; Actual arrow (tip) = " << vel / (rGrid->dx / RENDER_SIZE_CELL) << std::endl;

            //drawArrow(pos, vel, 0.8*cellSize);

            Vec2d point = vel;
            Vec2d base = pos;



            Vec2d w = point - base;
            double len = mag(w);

            if (len != 0) {
                w = w / (double) len; // normalize to build coordinate system

                // u = w + 90 
                // using rotation matrix  0  1
                //	                     -1  0
                Vec2d u = Vec2d(1 * w[1], -1 * w[0]);
                u = u / mag(u);

                // v = w - 90 (in fact v=-u)
                Vec2d v = Vec2d(-1 * w[1], 1 * w[0]);
                v = v / mag(v);


                double arrow_head_length = 0.05 * len;


                // arrow head points
                Vec2d arrow1, arrow2;
                arrow1 = point + arrow_head_length * (v - w);
                arrow2 = point + arrow_head_length * (u - w);


                std::cout << "base[0] = " << base[0] / (rGrid->dx / RENDER_SIZE_CELL) << ",  base[1] = " << base[1] / (rGrid->dx / RENDER_SIZE_CELL) << std::endl;
                std::cout << "point[0] = " << point[0] / (rGrid->dx / RENDER_SIZE_CELL) << ",  point[1] = " << point[1] / (rGrid->dx / RENDER_SIZE_CELL) << std::endl;
                std::cout << "point[0] = " << point[0] / (rGrid->dx / RENDER_SIZE_CELL) << ",  point[1] = " << point[1] / (rGrid->dx / RENDER_SIZE_CELL) << std::endl;
                std::cout << "arrow1[0] = " << arrow1[0] / (rGrid->dx / RENDER_SIZE_CELL) << ",  arrow1[1] = " << arrow1[1] / (rGrid->dx / RENDER_SIZE_CELL) << std::endl;
                std::cout << "point[0] = " << point[0] / (rGrid->dx / RENDER_SIZE_CELL) << ",  point[1] = " << point[1] / (rGrid->dx / RENDER_SIZE_CELL) << std::endl;
                std::cout << "arrow2[0] = " << arrow2[0] / (rGrid->dx / RENDER_SIZE_CELL) << ",  arrow2[1] = " << arrow2[1] / (rGrid->dx / RENDER_SIZE_CELL) << std::endl;
                std::cout << "--RENDERING--" << std::endl;
                std::cout << "rendering (base)= " << base / (rGrid->dx / RENDER_SIZE_CELL);
                std::cout << "; rendering (point)= " << point / (rGrid->dx / RENDER_SIZE_CELL) << std::endl;

            } else {
                std::cout << "--RENDERING CASE ZERO--" << std::endl;
                std::cout << "rendering (base)= " << base / (rGrid->dx / RENDER_SIZE_CELL);
                std::cout << "; rendering (point)= " << point / (rGrid->dx / RENDER_SIZE_CELL) << std::endl;
            }

        }
    return true;
}

/**
 * 
 */
void assignTestValues(Grid* rGrid) {

    // remember..
    // ni is the dimension in y axis.. the rows .. yes
    // nj is the dimension in x axis.. the columns.. hm
    // u has extra info in nj direction and v has extra info in ni direction


    for (int y = 0; y < rGrid->ni - 1; y++)
        for (int x = rGrid->nj / 2; x < rGrid->nj - 1; x++) {

            rGrid->u(x,y) = 30.0;
            rGrid->u(x + 1, y)=  30.0;
            rGrid->u(x, y + 1)= 30.0;
            rGrid->u(x + 1, y + 1) = 30.0;

            rGrid->v(x, y) = 50.0;
            rGrid->v(x + 1, y) = 50.0;
            rGrid->v(x + 1, y + 1) = 50.0;
            rGrid->v(x, y + 1) =50.0;

        }


    rGrid->pressure(0,0) = 30.0;
    rGrid->pressure(rGrid->ni-1, rGrid->nj-1) = 30;
    rGrid->pressure(rGrid->ni / 2, rGrid->nj / 2) = 30.0;


}

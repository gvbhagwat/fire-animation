/**
 * @file: RenderingInterface.cpp
 *
 * @date: 11-Jun-2012
 * @author: gaurav
 */

#include "Renderer.h"
#include "sim.h"
#include "MetaConfig.h"

//implicit inclusion of the Grid class from Grid.hpp

#include <cmath>
#include <GL/glut.h>
#include <GL/glu.h>
#include <iostream>

/**
 * 
 * @return cellSize of the grid
 */
int Renderer::getCellSize() {
    return cellSize;
}

/**
 * constructor
 */
Renderer::Renderer() {
    //  Auto-generated constructor stub
    if (META_LOG)
        std::cout << "--META--\tConstructor\tRenderer" << std::endl;

    cellSize = RENDER_SIZE_CELL;

    this->optionDrawGridCells = false;
    this->optionDrawFluidBody = true;

    this->optionDrawFlipParticles = false;
    this->optionHighlightFluidCells = false;
	this->optionHighlightFluidBoundaryCells = false;
    this->optionDrawPressure = false;
    this->optionDrawVelocity = false;
    this->optionDrawSimBoundary = false;
    this->optionDisplayDivergence = false;
    
    this->optionDiagnoseCustomDivergence_uComp = false;
    this->optionDiagnoseCustomDivergence_vComp = false;
    this->optionDiagnoseDivergence_uComp = false;
    this->optionDiagnoseDivergence_vComp = false;
}

/**
 * destructor
 */
Renderer::~Renderer() {
    //  Auto-generated destructor stub
    if (META_LOG)
        std::cout << "--META--\tDestructor\tRenderer" << std::endl;
}

/**
 *
 *
 */
void Renderer::drawSimulationEnitites(Grid* rGrid) {
    if (META_LOG)
        std::cout << "--META--\tFunctionCall\tRenderer::drawSimBoundary"
            << std::endl;

    for (unsigned int i = 0; i < rGrid->fireParticles.size(); i++) {
        Vec2d pos(rGrid->fireParticles[i]->pos / rGrid->dx);

        int x = int(pos[0]);
        int y = int(pos[1]);

        rGrid->marker(y, x) = 1; // means fluid cells
    }

    double red = 0.0, green = 0.0, blue = 0.0;

    for (int i = 0; i < rGrid->ni; i++) {
        for (int j = 0; j < rGrid->nj; j++) {

            if (rGrid->marker(i, j) == 2)
                red = green = blue = 0.75; //GRAY

            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glBegin(GL_POLYGON);
            glColor3d(red, green, blue);
            glVertex3d(j * (cellSize), (cellSize + cellSize * i), 0.0);
            glVertex3d(j * (cellSize), (0.0 + cellSize * i), 0.0);
            glVertex3d((j + 1)*(cellSize), (0.0 + cellSize * i), 0.0);
            glVertex3d((j + 1)*(cellSize), (cellSize + cellSize * i), 0.0);
            glEnd();

            red = green = blue = 0.0;
        }
    }

    this->drawMarkerParticles(rGrid);

    if (RENDER_LOG)
        std::cout << "--RENDER--\tdomain Boundary NOT drawn yet" << std::endl;
}

/**
 * Simply draws the grid according to the RENDER_SIZE cell size that specifies
 * the actual size of a cell for drawing. The granularity of grid is determined
 * by the a suitable variable from the Grid class.
 *
 * For current implementation, I have used grid->gridSize as the variable
 * for number of cells
 *
 * the drawing size of cell is different determined by the RENDER_SIZE
 *
 * @param grid, the object of the Grid class. For the current method
 * we only need the number of cells inside the grid which can be taken
 * from Grid class.
 *
 * TODO port it to OpengGL 4.0 or higher
 *
 * @see Grid Class
 */
void Renderer::drawGridCells(Grid* rGrid) {
    if (META_LOG)
        std::cout << "--META--\tFunctionCall\tRenderer::drawGridCells"
            << std::endl;

    int ni = rGrid->ni;
    int nj = rGrid->nj;

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++) {
            glBegin(GL_POLYGON);
            glColor3d(1.0, 1.0, 1.0);
            glVertex3d(j * (cellSize), (cellSize + cellSize * i), 0.0);
            glVertex3d(j * (cellSize), (0.0 + cellSize * i), 0.0);
            glVertex3d((j + 1)*(cellSize), (0.0 + cellSize * i), 0.0);
            glVertex3d((j + 1)*(cellSize), (cellSize + cellSize * i), 0.0);
            glEnd();
        }
    }

    if (META_LOG)
        std::cout << "--META--\tFunctionCall\tRenderer::drawGridCells"
            << std::endl;
}

/**
 * Simply draws the grid according to the RENDER_SIZE cell size that specifies
 * the actual size of a cell for drawing. The granularity of grid is determined
 * by the a suitable variable from the Grid class.
 *
 * For current implementation, I have used grid->gridSize as the variable
 * for number of cells
 *
 * the drawing size of cell is different determined by the RENDER_SIZE
 *
 * @param grid, the object of the Grid class. For the current method
 * we only need the number of cells inside the grid which can be taken
 * from Grid class.
 *
 * @todo port it to OpengGL 4.0 or higher
 *
 * @see Grid Class
 */
void Renderer::drawFluidBody(Grid* rGrid) {
    if (META_LOG)
        std::cout << "--META--\tFunctionCall\tRenderer::drawFluidBody" << std::endl;

    drawMarkerParticles(rGrid);

    if (RENDER_LOG)
        std::cout << "--RENDER--\tfluid body NOT drawn yet" << std::endl;
}

/**
 * @param grid
 * Simply draws the grid according to the RENDER_SIZE cell size that specifies
 * the actual size of a cell for drawing. The granularity of grid is determined
 * by the a suitable variable from the Grid class.
 *
 * For current implementation, I have used grid->gridSize as the variable
 * for number of cells
 *
 * the drawing size of cell is different determined by the RENDER_SIZE
 *
 * @param grid, the object of the Grid class. For the current method
 * we only need the number of cells inside the grid which can be taken
 * from Grid class.
 *
 * TODO port it to OpengGL 4.0 or higher
 *
 * @see Grid Class
 */
void Renderer::drawPressure(Grid* rGrid) {
    if (META_LOG)
        std::cout << "--META--\tFunctionCall\tRenderer::drawPressure" << std::endl;
    int ni = rGrid->ni;
    int nj = rGrid->nj;

    double maxPressure = 0.0;
    for (int i = 0; i < rGrid->ni; i++)
        for (int j = 0; j < rGrid->nj; j++)
            maxPressure = (maxPressure > rGrid->temp_pressure(i, j) ? maxPressure : rGrid->temp_pressure(i, j));

    //cout<<"max Pressure is = "<<maxPressure<<endl;

    double red = 0.0, green = 0.0, blue = 0.0;

    for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++) {


            if (rGrid->temp_pressure(i, j) > 0.0 && rGrid->temp_pressure(i, j) < maxPressure / 3.0) {
                red = 1.0;
                green = 0.0;
                blue = 0.0;
            } else if (rGrid->temp_pressure(i, j) >= maxPressure / 3.0 && rGrid->temp_pressure(i, j) < maxPressure * 2.0 / 3.0) {
                red = .0;
                green = 1.0;
                blue = 0.0;
            } else if (rGrid->temp_pressure(i, j) >= maxPressure * 2.0 / 3.0) {

                red = 0.0;
                green = 0.0;
                blue = 1.0;
            }

            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glBegin(GL_POLYGON);
            glColor3d(red, green, blue);
            glVertex3d(j * (cellSize), (cellSize + cellSize * i), 0.0);
            glVertex3d(j * (cellSize), (0.0 + cellSize * i), 0.0);
            glVertex3d((j + 1)*(cellSize), (0.0 + cellSize * i), 0.0);
            glVertex3d((j + 1)*(cellSize), (cellSize + cellSize * i), 0.0);
            glEnd();

            red = green = blue = 0.0;
        }
    }
    if (RENDER_LOG)
        std::cout << "--RENDER--\tpressure NOT drawn yet" << std::endl;
}


/**************************************************************************
 * THE divergence function and related diagnostics
 **************************************************************************
 * 1. displayDivergence
 * 2. displayCustomDivergence
 * 3. displayDiagnoseCustomDivergence_uComp(Grid* rGrid);
 * 4. displayDiagnoseCustomDivergence_vComp(Grid* rGrid);
 * 5. displayDiagnoseDivergence_vComp(Grid* rGrid);
 * 6. displayDiagnoseDivergence_vComp(Grid* rGrid);
 **************************************************************************/
void Renderer::displayDivergence(Grid* rGrid) {
    if (META_LOG)
        std::cout << "--META--\tFunctionCall\tRenderer::displayCustomDivergence" << std::endl;


    int ni = rGrid->ni;
    int nj = rGrid->nj;

    double red = 0.0, green = 0.0, blue = 0.0;

    for (int y = 0; y < ni; y++) {
        for (int x = 0; x < nj; x++) {
            if (rGrid->rhs(y, x) < 0.0) {
                red = 1.0;
                green = 1.0;
                blue = 0.0;
            }

            if (rGrid->rhs(y, x) > 0.0) {
                red = 1.0;
                green = 0.0;
                blue = 0.0;
            }
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glBegin(GL_POLYGON);
            glColor3d(red, green, blue);
            glVertex3d(x * (cellSize), (cellSize + cellSize * y), 0.0);
            glVertex3d(x * (cellSize), (0.0 + cellSize * y), 0.0);
            glVertex3d((x + 1)*(cellSize), (0.0 + cellSize * y), 0.0);
            glVertex3d((x + 1)*(cellSize), (cellSize + cellSize * y), 0.0);
            glEnd();

            red = green = blue = 0.0;
        }
    }
    
    if (RENDER_LOG)
        std::cout << "--RENDER--\trhs NOT drawn yet" << std::endl;


}

void Renderer::displayDiagnoseDivergence_uComp(Grid* rGrid){
        if (META_LOG)
        std::cout << "--META--\tFunctionCall\tRenderer::displayCustomDivergence" << std::endl;


    int ni = rGrid->ni;
    int nj = rGrid->nj;

    double red = 0.0, green = 0.0, blue = 0.0;

    for (int y = 0; y < ni; y++) {
        for (int x = 0; x < nj; x++) {
            if (rGrid->rhs_u_comp(y, x) < 0.0) {
                red = 0.0;
                green = 1.0;
                blue = 0.0;
            }

            if (rGrid->rhs_u_comp(y, x) > 0.0) {
                red = 0.0;
                green = 0.0;
                blue = 1.0;
            }
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glBegin(GL_POLYGON);
            glColor3d(red, green, blue);
            glVertex3d(x * (cellSize), (cellSize + cellSize * y), 0.0);
            glVertex3d(x * (cellSize), (0.0 + cellSize * y), 0.0);
            glVertex3d((x + 1)*(cellSize), (0.0 + cellSize * y), 0.0);
            glVertex3d((x + 1)*(cellSize), (cellSize + cellSize * y), 0.0);
            glEnd();

            red = green = blue = 0.0;
        }
    }
    
    if (RENDER_LOG)
        std::cout << "--RENDER--\trhs NOT drawn yet" << std::endl;

    
}

void Renderer::displayDiagnoseDivergence_vComp(Grid* rGrid){
        if (META_LOG)
        std::cout << "--META--\tFunctionCall\tRenderer::displayCustomDivergence" << std::endl;


    int ni = rGrid->ni;
    int nj = rGrid->nj;

    double red = 0.0, green = 0.0, blue = 0.0;

    for (int y = 0; y < ni; y++) {
        for (int x = 0; x < nj; x++) {
            if (rGrid->rhs_v_comp(y, x) < 0.0) {
                red = 0.0;
                green = 1.0;
                blue = 0.0;
            }

            if (rGrid->rhs_v_comp(y, x) > 0.0) {
                red = 0.0;
                green = 0.0;
                blue = 1.0;
            }
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glBegin(GL_POLYGON);
            glColor3d(red, green, blue);
            glVertex3d(x * (cellSize), (cellSize + cellSize * y), 0.0);
            glVertex3d(x * (cellSize), (0.0 + cellSize * y), 0.0);
            glVertex3d((x + 1)*(cellSize), (0.0 + cellSize * y), 0.0);
            glVertex3d((x + 1)*(cellSize), (cellSize + cellSize * y), 0.0);
            glEnd();

            red = green = blue = 0.0;
        }
    }
    
    if (RENDER_LOG)
        std::cout << "--RENDER--\trhs NOT drawn yet" << std::endl;    
}

void Renderer::displayCustomDivergence(Grid* rGrid) {
    if (META_LOG)
        std::cout << "--META--\tFunctionCall\tRenderer::displayCustomDivergence" << std::endl;

    int ni = rGrid->ni;
    int nj = rGrid->nj;

    double red = 0.0, green = 0.0, blue = 0.0;

    for (int y = 0; y < ni; y++)
        for (int x = 0; x < nj; x++) {
            if (rGrid->customRhs(y, x) < 0.0) {
                red = 1.0;
                green = 1.0;
                blue = 0.0;
            }
            if (rGrid->customRhs(y, x) > 0.0) {
                red = 1.0;
                green = 0.0;
                blue = 0.0;
            }

            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glBegin(GL_POLYGON);
            glColor3d(red, green, blue);
            glVertex3d(x * (cellSize), (cellSize + cellSize * y), 0.0);
            glVertex3d(x * (cellSize), (0.0 + cellSize * y), 0.0);
            glVertex3d((x + 1)*(cellSize), (0.0 + cellSize * y), 0.0);
            glVertex3d((x + 1)*(cellSize), (cellSize + cellSize * y), 0.0);
            glEnd();

            red = green = blue = 0.0;
        }


    if (RENDER_LOG)
        std::cout << "--RENDER--\trhs NOT drawn yet" << std::endl;



}

void Renderer::displayDiagnoseCustomDivergence_uComp(Grid* rGrid){
        if (META_LOG)
        std::cout << "--META--\tFunctionCall\tRenderer::displayCustomDivergence" << std::endl;


    int ni = rGrid->ni;
    int nj = rGrid->nj;

    double red = 0.0, green = 0.0, blue = 0.0;

    for (int y = 0; y < ni; y++) {
        for (int x = 0; x < nj; x++) {
            if (rGrid->diagnose_u_rhs(y, x) < 0.0) {
                red = 0.0;
                green = 1.0;
                blue = 0.0;
            }

            if (rGrid->diagnose_u_rhs(y, x) > 0.0) {
                red = 0.0;
                green = 0.0;
                blue = 1.0;
            }
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glBegin(GL_POLYGON);
            glColor3d(red, green, blue);
            glVertex3d(x * (cellSize), (cellSize + cellSize * y), 0.0);
            glVertex3d(x * (cellSize), (0.0 + cellSize * y), 0.0);
            glVertex3d((x + 1)*(cellSize), (0.0 + cellSize * y), 0.0);
            glVertex3d((x + 1)*(cellSize), (cellSize + cellSize * y), 0.0);
            glEnd();

            red = green = blue = 0.0;
        }
    }
}

void Renderer::displayDiagnoseCustomDivergence_vComp(Grid* rGrid){
        if (META_LOG)
        std::cout << "--META--\tFunctionCall\tRenderer::displayCustomDivergence" << std::endl;


    int ni = rGrid->ni;
    int nj = rGrid->nj;

    double red = 0.0, green = 0.0, blue = 0.0;

    for (int y = 0; y < ni; y++) {
        for (int x = 0; x < nj; x++) {
            if (rGrid->diagnose_v_rhs(y, x) < 0.0) {
                red = 0.0;
                green = 1.0;
                blue = 0.0;
            }

            if (rGrid->diagnose_v_rhs(y, x) > 0.0) {
                red = 0.0;
                green = 0.0;
                blue = 1.0;
            }
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glBegin(GL_POLYGON);
            glColor3d(red, green, blue);
            glVertex3d(x * (cellSize), (cellSize + cellSize * y), 0.0);
            glVertex3d(x * (cellSize), (0.0 + cellSize * y), 0.0);
            glVertex3d((x + 1)*(cellSize), (0.0 + cellSize * y), 0.0);
            glVertex3d((x + 1)*(cellSize), (cellSize + cellSize * y), 0.0);
            glEnd();

            red = green = blue = 0.0;
        }
    }
    
    if (RENDER_LOG)
        std::cout << "--RENDER--\trhs NOT drawn yet" << std::endl;    
}



/**************************************************************************
 * THE divergence function and related diagnostics
 **************************************************************************
 * 1. displayDivergence
 * 2. displayCustomDivergence
 * 3. displayDiagnoseCustomDivergence_uComp(Grid* rGrid);
 * 4. displayDiagnoseCustomDivergence_vComp(Grid* rGrid);
 * 5. displayDiagnoseDivergence_vComp(Grid* rGrid);
 * 6. displayDiagnoseDivergence_vComp(Grid* rGrid);
 **************************************************************************/


void Renderer::drawMarkerParticles(Grid* rGrid) {

    double scaleFactor = rGrid->dx / cellSize;
    
    for (unsigned int i = 0; i < rGrid->fireParticles.size(); i++) {
        //Vec2d pt = rGrid->markerParticles[i];
        
        //if (rGrid->fireParticles[i]->timeAlive < FLAME_HEIGHT){
        Vec2d pt = rGrid->fireParticles[i]->pos;
        pt = pt / scaleFactor;
        glPointSize(2);
        glBegin(GL_POINTS);
        if(rGrid->fireParticles[i]->timeAlive > 0.0 && rGrid->fireParticles[i]->timeAlive < 0.1)
              glColor3f(1.0,0.0,0.0);
        else if (rGrid->fireParticles[i]->timeAlive > 0.1 && rGrid->fireParticles[i]->timeAlive < 0.18)
              glColor3f(255.0/255.0, 140.0/255.0, 0.0);
        else
            glColor3f(1.0,1.0,0.0);
        glVertex2d(pt[0], pt[1]);
        glEnd();
        //}
    }
    

    glPopAttrib();
}

void Renderer::drawFlipParticles(Grid* rGrid) {

    double scaleFactor = rGrid->dx / cellSize;

    for (unsigned int i = 0; i < rGrid->flipParticles.size(); i++) {
        Vec2d pt = rGrid->flipParticles[i];
        pt = pt / scaleFactor;
        glPointSize(2);
        glBegin(GL_POINTS);
        glColor3f(1.0, 1.0, 0.0);
        glVertex2d(pt[0], pt[1]);
        glEnd();
    }

    glPopAttrib();


}


/**
 * Simply draws the grid according to the RENDER_SIZE cell size that specifies
 * the actual size of a cell for drawing. The granularity of grid is determined
 * by the a suitable variable from the Grid class.
 *
 * For current implementation, I have used grid->gridSize as the variable
 * for number of cells
 *
 * the drawing size of cell is different determined by the RENDER_SIZE
 *
 * @param grid, the object of the Grid class. For the current method
 * we only need the number of cells inside the grid which can be taken
 * from Grid class.
 *
 * TODO port it to OpengGL 4.0 or higher
 * TODO division is currently =4, need to develop a better color scheme
 *
 * @see Grid Class
 */
void Renderer::drawResultantVelocity(Grid* rGrid) {
    if (META_LOG)
        std::cout << "--META--\tFunctionCall\tRenderer::drawResultantVelocity"
            << std::endl;

    for (int y = 0; y < rGrid->nj; ++y) for (int x = 0; x < rGrid->ni; ++x) {

            Vec2d pos((x + 0.5) * rGrid->dx, (y + 0.5) * rGrid->dx);
            Vec2d vel = rGrid->getVelocity(pos);
            double colorCode[3];

            double magnitude = mag(vel);
            double maxMagnitude = mag(rGrid->maxVelocity);

            if (magnitude > maxMagnitude) {
                rGrid->maxVelocity = vel;
            }

            if (magnitude) // check the velocity is not Zero
            {
                vel = vel / mag(vel);
                vel = vel * rGrid->dx * 0.5;
            }
            vel = pos + vel;

            if (maxMagnitude > 1e-12 && magnitude < maxMagnitude / 3.0) {
                colorCode[0] = 1.0;
                colorCode[1] = 0.0;
                colorCode[2] = 0.0;
            } else if (maxMagnitude > 1e-12 &&
                    magnitude >= maxMagnitude / 3.0 &&
                    magnitude < maxMagnitude * (2.0) / 3.0) {
                colorCode[0] = 0.0;
                colorCode[1] = 1.0;
                colorCode[2] = 0.0;
            }
            if (maxMagnitude > 1e-12 &&
                    magnitude > maxMagnitude * (2.0) / 3.0) {
                colorCode[0] = 1.0;
                colorCode[1] = 0.0;
                colorCode[2] = 1.0;
            }

            double scaleFactor = rGrid->dx / cellSize;
            drawArrow(pos, vel, 0.00025 * cellSize, colorCode, scaleFactor);

            colorCode[0] = colorCode[1] = colorCode[2] = 1.0;

        }

    if (RENDER_LOG)
        std::cout << "--RENDER--\tvelocity vectors drawn" << std::endl;

}

/**
 * 
 * @param base
 * @param point
 * @param arrow_head_length
 */
void Renderer::drawArrow(Vec2d start, Vec2d end, double arrow_head_len, double colorCode[3], double scaleFactor) {
    if (META_LOG)
        std::cout << "--META--\tFunctionCall\tRenderer::drawArrow"
            << std::endl;

    Vec2d direction = end - start;

    Vec2d dir_norm = direction;

    if (mag(dir_norm) < 1e-14)
        return;

    normalize(dir_norm);
    Vec2d perp(dir_norm[1], -dir_norm[0]);

    Vec2d tip_left = end + arrow_head_len / (float) sqrt(2.0)*(-dir_norm + perp);
    Vec2d tip_right = end + arrow_head_len / (float) sqrt(2.0)*(-dir_norm - perp);



    start = start / scaleFactor;
    end = end / scaleFactor;
    tip_left = tip_left / scaleFactor;
    tip_right = tip_right / scaleFactor;


    glBegin(GL_LINES);
    glColor3dv(colorCode);
    glVertex2dv(start.v);
    glVertex2dv(end.v);
    glVertex2dv(end.v);
    glVertex2dv(tip_left.v);
    glVertex2dv(end.v);
    glVertex2dv(tip_right.v);
    glEnd();

    glPopAttrib();

}

/**
 * 
 * @param base
 * @param point
 * @param arrow_head_length
 */
void Renderer::drawArrow(Vec2d start, Vec2d end, double arrow_head_len, double scaleFactor) {
    if (META_LOG)
        std::cout << "--META--\tFunctionCall\tRenderer::drawArrow"
            << std::endl;

    Vec2d direction = end - start;

    Vec2d dir_norm = direction;

    //TODO Possibly automatically scale arrowhead length based on vector magnitude
    if (mag(dir_norm) < 1e-14)
        return;

    normalize(dir_norm);
    Vec2d perp(dir_norm[1], -dir_norm[0]);

    Vec2d tip_left = end + arrow_head_len / (float) sqrt(2.0)*(-dir_norm + perp);
    Vec2d tip_right = end + arrow_head_len / (float) sqrt(2.0)*(-dir_norm - perp);



    start = start / scaleFactor;
    end = end / scaleFactor;
    tip_left = tip_left / scaleFactor;
    tip_right = tip_right / scaleFactor;


    glBegin(GL_LINES);
    glColor3d(1.0, 0.0, 0.0);
    glVertex2dv(start.v);
    glVertex2dv(end.v);
    glVertex2dv(end.v);
    glVertex2dv(tip_left.v);
    glVertex2dv(end.v);
    glVertex2dv(tip_right.v);
    glEnd();

    glPopAttrib();

}

void Renderer::drawCircle(const Vec2f& centre, float rad, int segs) {
    glBegin(GL_POLYGON);
    for (int i = 0; i < segs; i++) {
        float cosine = rad * cos(i * 2 * M_PI / (float) (segs));
        float sine = rad * sin(i * 2 * M_PI / (float) (segs));
        glVertex2fv((Vec2f(cosine, sine) + centre).v);
    }
    glEnd();
}

void Renderer::highlightFluidCells(Grid* rGrid) {

    for (unsigned int i = 0; i < rGrid->fireParticles.size(); i++) {
        Vec2d pos(rGrid->fireParticles[i]->pos / rGrid->dx);

        int x = int(pos[0]);
        int y = int(pos[1]);

        rGrid->marker(y, x) = 1; // means fluid cells
    }

    double red = 0.0, green = 0.0, blue = 0.0;

    for (int i = 0; i < rGrid->ni; i++) {
        for (int j = 0; j < rGrid->nj; j++) {

            if (rGrid->marker(i, j) == 1) {
                red = green = blue = 1.0;
                // cout<<"found at i = "<<i<<" and j = "<<j<<endl;
            }

            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glBegin(GL_POLYGON);
            glColor3d(red, green, blue);
            glVertex3d(j * (cellSize), (cellSize + cellSize * i), 0.0);
            glVertex3d(j * (cellSize), (0.0 + cellSize * i), 0.0);
            glVertex3d((j + 1)*(cellSize), (0.0 + cellSize * i), 0.0);
            glVertex3d((j + 1)*(cellSize), (cellSize + cellSize * i), 0.0);
            glEnd();

            red = green = blue = 0.0;
        }
    }
    if (RENDER_LOG)
        std::cout << "--RENDER--\tpressure NOT drawn yet" << std::endl;

}

void Renderer::highlightFluidBoundaryCells(Grid* rGrid){

    double red = 0.0, green = 0.0, blue = 0.0;

    for (int i = 0; i < rGrid->ni; i++) {
        for (int j = 0; j < rGrid->nj; j++) {

            if (rGrid->boundary(i, j) == 1) {
                red = green = blue = 1.0;
                // cout<<"found at i = "<<i<<" and j = "<<j<<endl;
            }

            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glBegin(GL_POLYGON);
            glColor3d(red, green, blue);
            glVertex3d(j * (cellSize), (cellSize + cellSize * i), 0.0);
            glVertex3d(j * (cellSize), (0.0 + cellSize * i), 0.0);
            glVertex3d((j + 1)*(cellSize), (0.0 + cellSize * i), 0.0);
            glVertex3d((j + 1)*(cellSize), (cellSize + cellSize * i), 0.0);
            glEnd();

            red = green = blue = 0.0;
        }
    }
    if (RENDER_LOG)
        std::cout << "--RENDER--\tpressure NOT drawn yet" << std::endl;

}

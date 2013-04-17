/**
 * RenderingInterface.hpp
 *
 * @version: 1.0
 * @date: 11-Jun-2012
 * @author: gaurav
 */

#ifndef RENDERINGINTERFACE_HPP_
#define RENDERINGINTERFACE_HPP_

#include "Grid.h"
#include "Vec.h"

/**
 * Rendering class handles the OpenGL rendering code modeled as an interface.
 *
 * Current Implementation only handles opengl 2.1, will port to
 * OpenGL version 4.0 soon.
 */
class Renderer {
//protected:

    //int cellSize;

public:

    int cellSize;
    bool optionDrawVelocity;
    bool optionDrawGridCells;
    bool optionDrawSimBoundary;
    bool optionDrawPressure;
    bool optionDrawFluidBody;
    bool optionDrawFlipParticles;
    bool optionHighlightFluidCells;
    bool optionDisplayDivergence;
    bool optionDisplayCustomDivergence;
    
    // diagnostic purposes
    bool optionDiagnoseCustomDivergence_uComp;
    bool optionDiagnoseCustomDivergence_vComp;
    
    bool optionDiagnoseDivergence_uComp;
    bool optionDiagnoseDivergence_vComp;
    
    Renderer();
    virtual
    ~Renderer();

    int getCellSize();
    
    void drawGridCells(Grid* rGrid);

    void drawSimulationEnitites(Grid* rGrid);

    void drawResultantVelocity(Grid* rGrid);

    void drawFluidBody(Grid* rGrid);

    void drawPressure(Grid* rGrid);
    
    void highlightFluidCells(Grid* rGrid);
    
    void displayDivergence(Grid* grid);
    
    void displayDiagnoseCustomDivergence_uComp(Grid* rGrid);
    void displayDiagnoseCustomDivergence_vComp(Grid* rGrid);
    
    void displayDiagnoseDivergence_uComp(Grid* rGrid);
    void displayDiagnoseDivergence_vComp(Grid* rGrid);
    
    void displayCustomDivergence(Grid* grid);
    
    //TODO only made for this grid
    
    void drawMarkerParticles(Grid* grid);
    void drawFlipParticles(Grid* grid);
    
    
    // helper functions
    void drawArrow(Vec2d base, Vec2d point, double arrow_head_length, double scaleFactor);
    void drawArrow(Vec2d base, Vec2d point, double arrow_head_length, double colorCode[3], double scaleFactor);
    void drawCircle(const Vec2f& centre, float rad, int segs);
    

    
};

#endif /* RENDERINGINTERFACE_HPP_ */


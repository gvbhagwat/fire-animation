/**
 *  FLUID SOLVER BY : GAURAV BHAGWAT
 *
 */

// Compulsory Headers
#include <GL/glut.h>

#include <iostream>
#include <fstream>
#include <cmath>

// Solver Headers
#include "Renderer.h"
#include "Grid.h"
#include "GridTest.h"
#include "FluidSim.h"
#include "Dimensions.h"
#include "MetaConfig.h"

using namespace std;


// Frame Capturing Variables and declarations

void capture_frame(unsigned int framenum);

int SCREEN_WIDTH= GRID_NJ * PHYSICAL_WIDTH ;
int SCREEN_HEIGHT= GRID_NI * PHYSICAL_HEIGHT;


bool recording=false;
unsigned int framenum=0;
unsigned char *pRGB;

// Global Variables SIMULATION related
Grid* grid;
Renderer* renderObj;
FluidSim* sim;

FluidBodyInitialization choice = DAM_BREAK; // defined in FluidSim.h

double timestep = 0.01;

// OpenGL call back declarations
void init();
void display();
void reshape(int, int);
void keyboard(unsigned char, int, int);

// Testing Functions
void initGridTest();

/**
 *
 */
void init() {

    // create the grid
    grid = new Grid(GRID_NI, GRID_NJ, PHYSICAL_WIDTH, PHYSICAL_HEIGHT);

    // create the simulator
    sim = new FluidSim(*grid);

    // initialize the simulator, (also initializes the grid)
    sim->initialize(START_TIME, END_TIME);

    // initialize the solidBoundariues
    sim->initializeSolidBoundaries();

    // initialize the fluidBody
    sim->initializeFluidBody(2);

    // MOCK TEST
    // sim->advance(grid,timestep);

    // all that remains is to start the simulator
    // done by sim>advance() 
    // but its called by timerCallback

        // create the OpenGL renderer
    renderObj = new Renderer();

    //OpenGL setup
    //glShadeModel(GL_FLAT);
}

/**
 *
 *
 */
void cleanup() {
    //delete sim;
    delete renderObj;
    delete grid;
}

/**
 * 
 * @param junk
 */
void timer(int junk) {

    // advances the simulation by "timestep"
    sim->advance(timestep);

    glutPostRedisplay();
    glutTimerFunc(timestep, timer, 0);
}

/**
 *
 */
void display() {

    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
    glClearColor(0.0, 0.0, 0.0, 0.0);

    if (renderObj->optionDrawFluidBody) {
        renderObj->drawSimulationEnitites(grid);
    }

    if (renderObj->optionDrawVelocity) {
        renderObj->drawResultantVelocity(grid);
    }

    if (renderObj->optionDrawPressure) {
        renderObj->drawPressure(grid);
    }

    if (renderObj->optionHighlightFluidCells) {
        renderObj->highlightFluidCells(grid);
    }

    if (renderObj->optionDisplayDivergence) {
        renderObj->displayDivergence(grid);
    }

    if (renderObj->optionDisplayCustomDivergence) {
        renderObj->displayCustomDivergence(grid);
    }
    
    if(renderObj->optionDiagnoseDivergence_vComp){
        renderObj->displayDiagnoseDivergence_vComp(grid);
    }

    if(renderObj->optionDiagnoseDivergence_uComp){
        renderObj->displayDiagnoseDivergence_uComp(grid);
    }
    
    if(renderObj->optionDiagnoseCustomDivergence_vComp){
        renderObj->displayDiagnoseCustomDivergence_vComp(grid);
    }

    if(renderObj->optionDiagnoseCustomDivergence_uComp){
        renderObj->displayDiagnoseCustomDivergence_uComp(grid);
    }
    
    if (renderObj->optionDrawGridCells) {
        renderObj->drawGridCells(grid);
    }

    if(recording)
        capture_frame(framenum++);
    
    glutSwapBuffers();

}

/**
 *
 */
void positionCamera() {
    // Define the viewing matrix.
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslated(10.0, 10.0, 0.0);
}

/**
 *
 * @param w
 * @param h
 */
void reshape(int w, int h) {
    if (h == 0)
        h = 1;
    glViewport(0, 0, (GLsizei) w, (GLsizei) h);

    SCREEN_WIDTH = w;
    SCREEN_HEIGHT = h;
    
    // Define the projection matrix.
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, (GLdouble) w, 0.0, (GLdouble) h);

    positionCamera();
}

/**
 *
 * @param key
 * @param x
 * @param y
 */
void keyboard(unsigned char key, int x, int y) {
    switch (key) {
        case 'v':
            renderObj->optionDrawVelocity = true;
            break;
        case 'V':
            renderObj->optionDrawVelocity = false;
            break;
        case 'f':
            renderObj->optionDrawFluidBody = true;
            break;
        case 'F':
            renderObj->optionDrawFluidBody = false;
            break;
        case 'q':
            renderObj->optionDrawFlipParticles = true;
            break;
        case 'Q':
            renderObj->optionDrawFlipParticles = false;
            break;
        case 'h':
            renderObj->optionHighlightFluidCells = true;
            break;
        case 'H':
            renderObj->optionHighlightFluidCells = false;
            break;
        case 'p':
            renderObj->optionDrawPressure = true;
            break;
        case 'P':
            renderObj->optionDrawPressure = false;
            break;
        case 'g':
            renderObj->optionDrawGridCells = true;
            break;
        case 'G':
            renderObj->optionDrawGridCells = false;
            break;
        case 'c':
            renderObj->optionDisplayCustomDivergence = true;
            break;
        case 'C':
            renderObj->optionDisplayCustomDivergence= false;
            break;            
        case 'r':
            renderObj->optionDisplayDivergence = true;
            break;
        case 'R':
            renderObj->optionDisplayDivergence = false;
            break;
        case 'a':
            renderObj->optionDiagnoseDivergence_uComp = true;
            break;
        case 'A':
            renderObj->optionDiagnoseDivergence_uComp = false;
            break;            
        case 'z':
            renderObj->optionDiagnoseDivergence_vComp = true;
            break;
        case 'Z':
            renderObj->optionDiagnoseDivergence_vComp = false;
            break;                        
        case 's':
            renderObj->optionDiagnoseCustomDivergence_uComp = true;
            break;
        case 'S':
            renderObj->optionDiagnoseCustomDivergence_uComp = false;
            break;                                    
        case 'x':
            renderObj->optionDiagnoseCustomDivergence_vComp = true;
            break;
        case 'X':
            renderObj->optionDiagnoseCustomDivergence_vComp = false;
            break;                                                
        case 27:
            if (FULL_SCREEN) {
                glutLeaveGameMode();
            }
            cleanup();
            exit(0);
            break;
    }

    glutPostRedisplay();
    //glFlush();
}

/**
 *
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char* argv[]) {

    // Create all the objects and define the fluid
    init();

    if (CLASS_TESTING_MODE == true)
        initGridTest();
    else {

        //initGridTest();

        // Initialize the display window.
        glutInit(&argc, argv);
        glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);

        //defining the coordinate axes equal to the number_of_cells x size of each cell
        int width = grid->ni * renderObj->getCellSize();
        int height = grid->nj * renderObj->getCellSize();

        //The window
        glutInitWindowSize((width + 20), (height + 20)); //adding margin 20
        glutInitWindowPosition(width / 2, height / 2);
        glutCreateWindow("Fluid Solver - Gaurav Bhagwat");

        glutTimerFunc(1000, timer, 0);

        // Initialize the event handlers.
        glutDisplayFunc(display);
        glutReshapeFunc(reshape);
        glutKeyboardFunc(keyboard);

        // Start the program.
        glutMainLoop();


        exit(EXIT_SUCCESS);
    }

    return 0;
}

void initGridTest() {

    //checkInitialize(grid);
    //assignTestValues(grid);
    //checkGetVelocity(grid);
    //checkRenderingCoordinates(grid);

    cleanup();
}

void capture_frame(unsigned int framenum) {
    
    // int SCREEN_WIDTH = PHYSICAL_WIDTH * GRID_NJ +20;
    // int SCREEN_HEIGHT = PHYSICAL_HEIGHT * GRID_NI + 20;
    //global pointer float *pRGB
    pRGB = new unsigned char [3 * (SCREEN_WIDTH + 1) * (SCREEN_HEIGHT + 1) ];


    // set the framebuffer to read
    //default for double buffered
    glReadBuffer(GL_BACK);

    glPixelStoref(GL_PACK_ALIGNMENT, 1); //for word allignment

    glReadPixels(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT, GL_RGB, GL_UNSIGNED_BYTE, pRGB);
    char filename[200];
    sprintf(filename, "simframes/frame_%04d.ppm", framenum);
    std::ofstream out(filename, std::ios::out);
    out << "P6" << std::endl;
    out << SCREEN_WIDTH << " " << SCREEN_HEIGHT << " 255" << std::endl;
    out.write(reinterpret_cast<char const *> (pRGB), (3 * (SCREEN_WIDTH + 1) * (SCREEN_HEIGHT + 1)) * sizeof (int));
    out.close();

    //function to store pRGB in a file named count
    delete pRGB;
}


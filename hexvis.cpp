
#ifndef FLT
# error "Please define FLT when compiling (hint: See CMakeLists.txt)"
#endif


/*!
 * Include the calculation class.
 */
#include "calc.h"

/*!
 * Standard library includes here
 */
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <string>
#include <sstream>
#include <limits>
#include <chrono>
/*
 * Morphologica includes here
 */
# include <morph/Visual.h>
# include <morph/HexGridVisual.h>
# include <morph/ColourMap.h>
# include <morph/VisualDataModel.h>
# include <morph/Scale.h>
# include <morph/vec.h>
/*
 * Other useful morph includes.
 */
#include <morph/tools.h>
#include <morph/Config.h> //json read-writer

/*
main() is in crontrol of the simulation. Parameters are controlled via JSON,
argc and argv are used to pass in params.json at execution.*/


/*Lots of explict type definition is used in this code. e.g. 100UL and 0.3f*/

int main(int argc, char **argv){
    
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " /path/to/params.json" << std::endl;
        return 1;
    }
    std::string paramsfile (argv[1]);

    /*
     * Set up morph::Config JSON reader/writer for reading the parameters
     * Return error if this fails.
     */
    morph::Config conf(paramsfile);
    if (!conf.ready) {
        return 1;
    }
    

    //getting simulation-wide parameters from JSON
    // usage here: conf.get("find a value", if not found use this value:)
    const unsigned int steps = conf.getUInt ("steps", 1000UL);
    if (steps == 0) {
        std::cerr << "Not much point simulating 0 steps! Exiting." << std::endl;
        return 1;}

    const FLT dt = static_cast<FLT>(conf.getDouble("dt", 0.00001));

    //returning No of steps.
    std::cout << "steps to simulate: " << steps << std::endl;

/*
 *
 * Visualisation stuff.
 *
 */

// Parameters for plotting:
const unsigned int plotevery = conf.getUInt ("plotevery", 10);

// Window width and height
const unsigned int win_width = conf.getUInt ("win_width", 1025UL);
unsigned int win_height_default = static_cast<unsigned int>(0.8824f * (float)win_width);
const unsigned int win_height = conf.getUInt ("win_height", win_height_default);

// Set up the morph::Visual object:
morph::Visual v1 (win_width, win_height, "Hex Diffusion");

// Set a dark blue background, RGB:
v1.bgcolour = {conf.getFloat("bgR", 0.2f), conf.getFloat("bgG", 0.2f), conf.getFloat("bgB", 0.2f), 1.0f};
// setting clipping planes & FOV.
v1.zNear = 0.001;
v1.zFar = 10000;
v1.fov = 45;

// unlock the scene (you cna pan around the visual)
v1.sceneLocked = conf.getBool ("sceneLocked", false);

// You can set the scene x/y/z offsets.
v1.setZDefault (conf.getFloat ("z_default", 0.0f));
v1.setSceneTransXY (conf.getFloat ("x_default", 0.0f),
                    conf.getFloat ("y_default", 0.0f));

// scroll in-out speed.
v1.scenetrans_stepsize = 0.5;

// setting up the render clock.
std::chrono::steady_clock::time_point lastrender = std::chrono::steady_clock::now();



/*
*
*
* Instantiating & setting up the model object class found in calc.h
*
*
*/
    D_calc<FLT> D;

    D.svgpath = conf.getString ("svgpath", "");
    D.ellipse_a = conf.getDouble ("ellipse_a", 0.6);
    D.ellipse_b = conf.getDouble ("ellipse_b", 0.6);


    // Control the size of the hexes, and therefore the number of hexes in the grid
    D.hextohex_d = conf.getFloat ("hextohex_d", 0.01f);
    D.hexspan = conf.getFloat ("hexspan", 4.0f);
    //hexspan/hextohex = No of hexes.

    // Boundary fall-off distance, this is will ensure the system rolls off to zero at the boundary.
    //This is essentially neutrons being absorbed by the casing of the reactor.
    D.boundaryFalloffDist = conf.getFloat("boundaryFalloffDist", 0.1f);

    // After setting the first few features, we can call the allocate function to set
    D.allocate();

    // After allocate(), we can set up parameters:
    D.set_dt(dt);
    D.k1 = conf.getDouble ("k1", 1);
    D.k2 = conf.getDouble ("k2", 1);
    D.k3 = conf.getDouble ("k3", 1);
    D.k4 = conf.getDouble ("k4", 1);
    D.D_phi = conf.getDouble ("D_phi", 0.1);

    // Now parameters are set, call init(), which in this example simply initializes Phi with some noise.
    D.init();

/*
* This is the end of model setup.
*/



// Spatial offset, for positioning of HexGridVisuals
morph::vec<float> spatOff;
float xzero = 0.0f;

xzero -= 0.5 * D.hg->width();
spatOff = { 0.0, 0.0, 0.0 };


morph::ColourMapType cmt = morph::ColourMap<FLT>::strToColourMapType (conf.getString ("colourmap", "Jet"));

// Create a new HexGridVisual hgv1 then set its parameters:
auto hgv1 = std::make_unique<morph::HexGridVisual<FLT>> (D.hg, spatOff);
v1.bindmodel (hgv1);
hgv1->setScalarData (&D.phi);
// Z position scaling - how hilly/bumpy the visual will be.
hgv1->zScale.setParams (0.5f, 0.0f);
// The second is the colour scaling. Set this to autoscale.
hgv1->colourScale.do_autoscale = true;
hgv1->cm.setType (cmt);

// "finalize" is required before adding the HexGridVisual to the morph::Visual.
hgv1->finalize();
auto hgv1p = v1.addVisualModel (hgv1);


    // Start the loop
    bool finished = false;
    while (finished == false) {
        // Step the model
        D.step();

        if((D.stepCount % 100) == 0){
        std::cout << "My data range is " << D.phi.range() << std::endl;
        };

        if ((D.stepCount % plotevery) == 0) {
            // These two lines update the data for the two hex grids. That leads to
            // the CPU recomputing the OpenGL vertices for the visualizations.
            hgv1p->updateData (&(D.phi));
            hgv1p->clearAutoscaleColour();
        }
        // rendering the graphics. After each simulation step, check if enough time
        // has elapsed for it to be necessary to call v1.render().
        std::chrono::steady_clock::duration sincerender = std::chrono::steady_clock::now() - lastrender;
        if (std::chrono::duration_cast<std::chrono::milliseconds>(sincerender).count() > 17) { // 17 is about 60 Hz
            glfwPollEvents();
            v1.render();
            lastrender = std::chrono::steady_clock::now();
        }

        if (D.stepCount > steps) {
            finished = true;
        }
    }


    std::cout << "Ctrl-c or press x in graphics window to exit.\n";
    v1.keepOpen();

return 0;
}

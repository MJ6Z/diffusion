    //using FLT defined from Cmakelists so I can change the size of my floating point numbers as desired.
#ifndef FLT
# error "Please define FLT when compiling in CMakeLists.txt"
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
# include <mplot/Visual.h>
# include <mplot/HexGridVisual.h>
# include <mplot/HexGrid.h>
# include <mplot/ColourMap.h>
# include <mplot/VisualDataModel.h>
# include <sm/scale>
# include <sm/vec>
/*
* Other useful morph includes.
*/
#include <mplot/tools.h>
#include <sm/config> //json read-writer

/*
* main() is in control of the simulation. Parameters are controlled via JSON,
* argc and argv are used to pass in params.json at execution.*/


/*Lots of explict type definition is used in this code. e.g. 100UL and 0.3f*/

int main(int argc, char **argv){

    /*
    * Arguments
    */

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " /path/to/params.json" << std::endl;
        return 1;

    }
    std::string paramsfile (argv[1]);

    //Set up morph::Config JSON reader/writer for reading the parameters
    //Return error if this fails.

    morph::Config conf(paramsfile);
    if (!conf.ready) {
        return 1;

    }
    //get bool debug and option params for use later.
    //usage here: conf.get("find a value", if not found use this value:)
    bool debug = conf.getBool("debug",false);
    bool range_out = conf.getBool("range_out",false);
    bool doAutoScale = conf.getBool("doColourAutoscale",true);


    //get coordinate r,g,b positions for fuel rods
    morph::vvec<int> coolant_r = conf.getvvec<int>("coolant_r_locations");
    morph::vvec<int> coolant_g = conf.getvvec<int>("coolant_g_locations");
    morph::vvec<int> coolant_b = conf.getvvec<int>("coolant_b_locations");

    //check there are no incomplete coordinates.
    if(coolant_r.size() != coolant_g.size() || coolant_g.size() != coolant_b.size()){
        std::cerr<<"Ensure coolant positions are complete in paramsfile" <<std::endl;
        return 1;
    }
    //get coordinate r,g,b positions for contorl rods
    morph::vvec<int> control_r = conf.getvvec<int>("control_r_locations");
    morph::vvec<int> control_g = conf.getvvec<int>("control_g_locations");
    morph::vvec<int> control_b = conf.getvvec<int>("control_b_locations");

    //check there are no incomplete coordinates.
    if(control_r.size() != control_g.size() || control_g.size() != control_b.size()){
        std::cerr<<"Ensure control positions are complete in paramsfile" <<std::endl;
        return 1;
    }

    //get coordinate r,g,b positions for fuel rods.
    morph::vvec<int> fuel_r = conf.getvvec<int>("fuel_r_locations");
    morph::vvec<int> fuel_g = conf.getvvec<int>("fuel_g_locations");
    morph::vvec<int> fuel_b = conf.getvvec<int>("fuel_b_locations");

    //check there are no incomplete coordinates.
    if(fuel_r.size() != fuel_g.size() || fuel_g.size() != fuel_b.size()){
        std::cerr<<"Ensure fuel positions are complete in paramsfile" <<std::endl;
        return 1;
    }
    //get coordinate r,g,b positions for sources.
    morph::vvec<int> source_r = conf.getvvec<int>("source_r_locations");
    morph::vvec<int> source_g = conf.getvvec<int>("source_g_locations");
    morph::vvec<int> source_b = conf.getvvec<int>("source_b_locations");

    //check there are no incomplete coordinates.
    if(source_r.size() != source_g.size() || source_g.size() != source_b.size()){
        std::cerr<<"Ensure source positions are complete in paramsfile" <<std::endl;
        return 1;
    }


    //getting simulation-wide parameters from JSON

    //usage here: conf.get("find a value", if not found use this value:)
    const unsigned int steps = conf.getUInt ("steps", 1000UL);
    if (steps == 0) {
        std::cerr << "Not much point simulating 0 steps! Exiting." << std::endl;
        return 1;

    }
    //get the dt timestep value.
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
    v1.setSceneTransZ (conf.getFloat ("z_default", 0.0f));
    v1.setSceneTransXY (conf.getFloat ("x_default", 0.0f),
                        conf.getFloat ("y_default", 0.0f));

    // scroll in-out speed.
    v1.scenetrans_stepsize = 0.5;

    // setting up the render clock.
    std::chrono::steady_clock::time_point lastrender = std::chrono::steady_clock::now();



    /*
    *
    *
    * Instantiating & setting up the D object from the class D_calc found in calc.h
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

    // After setting the first few features, we can call the allocate function to set
    D.allocate();

    // After allocate(), we can set up parameters:

    //timestep
    D.set_dt(dt);

    //Runge Kutta parameters
    D.k1 = conf.getDouble ("k1", 1);
    D.k2 = conf.getDouble ("k2", 1);
    D.k3 = conf.getDouble ("k3", 1);
    D.k4 = conf.getDouble ("k4", 1);
    //Diffusion coefficients
    D.D_Fflux = conf.getDouble ("D_Fflux", 0.1);
    D.D_THflux = conf.getDouble ("D_THflux", 0.1);
    D.D_T = conf.getDouble ("D_T", 0.1);

    //Values on the hexgrid.
    //and doNoise being a option to create a random set of values to map to Fflux.
    D.doTemperatureNoise = conf.getBool ("doTemperatureNoise",false);
    D.doFluxNoise = conf.getBool ("doFluxNoise",false);
    D.noiseMinHeight = conf.getFloat("noiseMinHeight", 0);
    D.noiseMaxHeight = conf.getFloat("noiseMaxHeight", 0);


    D.sourceNeutrons = conf.getBool("sourceNeutrons",false);
    D.source_strength = conf.getDouble ("source_strength", 1);

    //the temperature at which the reactor will undergo a metltdown.
    D.temperature_meltdown_value = conf.getFloat("temperature_meltdown_value",1);


    //calculation coefficients
    D.dTdt_THflux_coeff = conf.getDouble ("dTdt_THflux_coeff", 0.1);
    D.dTdt_Fflux_coeff = conf.getDouble ("dTdt_Fflux_coeff", 0.1);
    D.neutrons_per_fission = conf.getDouble ("neutrons_per_fission", 0.1);
    D.fuel_strength = conf.getDouble ("fuel_strength", 0.1);



    D.moderation_strength = conf.getDouble ("moderation_strength", 0.1);
    if(D.moderation_strength>1){
        std::cerr<<"Cannot moderate more than 100 percent of neutrons - Exitting." << std::endl;
        return 1;
    };


    D.absorbtion_strength = conf.getDouble ("absorbtion_strength", 0.1);
    if(D.moderation_strength>1){
        std::cerr<<"Cannot absorb more than 100 percent of neutrons, check your parameters - Exitting." << std::endl;
        return 1;
    };


    // Passing in positions. See comments.
    {
        D.coolant_positions.resize(coolant_b.size());    //can use coolant_b.size as r,g,b have been validated to eb equal in size.
        for(unsigned int i=0; i<coolant_b.size(); i++){ //repacking data from individual coolant indexes into vec int, 3.
            D.coolant_positions[i]={coolant_r[i],coolant_g[i],coolant_b[i]};
        }


        D.control_positions.resize(control_b.size());    //can use control_b.size as r,g,b have been validated to eb equal in size.
        for(unsigned int i=0; i<control_b.size(); i++){ //repacking data from individual control indexes into vec int, 3.
            D.control_positions[i]={control_r[i],control_g[i],control_b[i]};
        }


        D.fuel_positions.resize(fuel_b.size());    //can use fuel_b.size as r,g,b have been validated to eb equal in size.
        for(unsigned int i=0; i<fuel_b.size(); i++){ //repacking data from individual fuel indexes into vec int, 3.
            D.fuel_positions[i]={fuel_r[i],fuel_g[i],fuel_b[i]};
        }

        D.source_positions.resize(source_b.size());    //can use source_b.size as r,g,b have been validated to eb equal in size.
        for(unsigned int i=0; i<source_b.size(); i++){ //repacking data from individual source indexes into vec int, 3.
            D.source_positions[i]={source_r[i],source_g[i],source_b[i]};
        }
    }


    // Now parameters are set, call init().
    D.init();

    /*
    * This is the end of model setup.
    */


    // Before starting the simulation, create the HexGridVisuals.

// scale params.
    float ZScaleMin;
    float ZScaleMax;
    ZScaleMin = conf.getFloat("ZScaleMin",0);
    ZScaleMax = conf.getFloat("ZScaleMax",0);

// Spatial offsets, for positioning of the HexGridVisuals.
    morph::vec<float> spatOff;
    float xzero = 0.0f;
    float yzero = 0.0f;
    xzero = D.hg->width();
    yzero = D.hg->width();

//Set the colourmaptype as various cmt_ varaibles. Note the use of strToColourMapType, a very useful function.
    morph::ColourMapType cmt_Fflux = morph::ColourMap<FLT>::strToColourMapType (conf.getString ("colourmap_Fflux", "Jet"));
    morph::ColourMapType cmt_THflux = morph::ColourMap<FLT>::strToColourMapType (conf.getString ("colourmap_THflux", "Jet"));
    morph::ColourMapType cmt_T = morph::ColourMap<FLT>::strToColourMapType (conf.getString ("colourmap_T", "Jet"));
    morph::ColourMapType cmt_celltype = morph::ColourMap<FLT>::strToColourMapType (conf.getString ("celltype_colourmap", "Jet"));

// Create a new HexGridVisual then set its parameters (zScale, colourScale, etc.
// this one is for Fflux.
    spatOff = { -0.9f*xzero, 0.0f, 0.0f };
    auto hgv1 = std::make_unique<morph::HexGridVisual<FLT>> (D.hg.get(), spatOff);
    v1.bindmodel (hgv1);
    hgv1->setScalarData (&D.Fflux);
    hgv1->zScale.setParams (ZScaleMin, ZScaleMax);

    // colour properties & labelling.
    hgv1->colourScale.do_autoscale = true;
    hgv1->cm.setType(cmt_Fflux);

    if(debug)
    {
        hgv1->addLabel("hgv1 binded to hgvp1 data=Fflux", { -0.2f, D.ellipse_b*-1.4f, 0.01f }, morph::TextFeatures(0.1f, morph::colour::white));
    }else{
        hgv1->addLabel("Fast neutron flux", { -0.2f, D.ellipse_b*-1.4f, 0.01f },
        morph::TextFeatures(0.1f, morph::colour::white));
    }

    // "finalize" is required before adding the HexGridVisual to the morph::Visual.
    hgv1->finalize();
    auto hgv1p = v1.addVisualModel (hgv1);

//temperature visual.
    spatOff = { xzero, 0.0f, 0.0f };
    auto hgv2 = std::make_unique<morph::HexGridVisual<FLT>> (D.hg.get(), spatOff);
    v1.bindmodel (hgv2);
    hgv2->setScalarData (&D.T);
    hgv2->zScale.setParams (ZScaleMin, ZScaleMax);
    hgv2->colourScale.do_autoscale = true;
    hgv2->cm.setType (cmt_T);
    if(debug){
        hgv2->addLabel ("hgv2 binded to hgvp2, data=D.T", { -0.2f, D.ellipse_b*-1.4f, 0.01f },
                        morph::TextFeatures(0.1f, morph::colour::white));
    }else{
        hgv2->addLabel ("Temperature scaled from max to min", { -0.2f, D.ellipse_b*-1.4f, 0.01f },
                        morph::TextFeatures(0.1f, morph::colour::white));
    }
    hgv2->finalize();
    auto hgv2p = v1.addVisualModel (hgv2);

//thermal flux visual.
    spatOff = {-0.9f*xzero, 1.3f*yzero, 0.0f };
    auto hgv3 = std::make_unique<morph::HexGridVisual<FLT>> (D.hg.get(), spatOff);
    v1.bindmodel (hgv3);
    hgv3->setScalarData (&D.THflux);
    hgv3->zScale.setParams (ZScaleMin, ZScaleMax);
    hgv3->colourScale.do_autoscale = true;
    hgv3->cm.setType (cmt_THflux);
    if(debug){
        hgv3->addLabel ("hgv3 binded to hgvp3, data=D.THflux", { -0.2f, D.ellipse_b*-1.4f, 0.01f },
                        morph::TextFeatures(0.1f, morph::colour::white));
    }else{
        hgv3->addLabel ("Thermal neuton flux", { -0.2f, D.ellipse_b*-1.4f, 0.01f },
                        morph::TextFeatures(0.1f, morph::colour::white));
    }
    hgv3->finalize();
    auto hgv3p = v1.addVisualModel(hgv3);


//type indicator visual.
    spatOff = {0.9f*xzero, 1.3f*yzero, 0.0f };
    auto hgv5 = std::make_unique<morph::HexGridVisual<FLT>> (D.hg.get() , spatOff);
    v1.bindmodel (hgv5);
    hgv5->setScalarData (&D.show_celltype);
    hgv5->zScale.setParams (ZScaleMin, ZScaleMax);
    hgv5->colourScale.do_autoscale = true;
    hgv5->cm.setType (cmt_celltype);
    hgv5->addLabel ("CONTROL(light blue), COOLANT(dark blue), FUEL (red)", { -0.2f, D.ellipse_b*-1.4f, 0.01f },
                    morph::TextFeatures(0.1f, morph::colour::white));
    hgv5->finalize();
    auto hgv5p = v1.addVisualModel(hgv5);





    // Start the loop
    bool finished = false;
    //while the loop isn't finished and ready to finish is false. ctrl+q triggers readytofinish to be true.
    while (finished == false && v1.readyToFinish == false) {

        if(D.stepCount == 1){
            hgv5p->updateData (&(D.show_celltype));
            hgv5p->clearAutoscaleColour();
        }

        if(range_out){ //if the rnge debug variable is true, spit out some values.
            if((D.stepCount % 10000) == 0 || D.stepCount == 0){
                std::cout << "Fflux.range = " << D.Fflux.range() << std::endl;
                std::cout << "THflux.range = " << D.THflux.range() << std::endl;
                std::cout << "t.range = " << D.T.range() << std::endl;


            }
        }

        // Step the model
        D.step();

        if (D.stepCount > steps) { //check if the program has finished, willl run last step.
            finished = true;
        }

        if ((D.stepCount % plotevery) == 0) {

            //updades data in the visuals only when the visual is about to be replotted to save resources.
            hgv1p->updateData (&(D.Fflux));
            hgv2p->updateData (&(D.T));
            hgv3p->updateData (&(D.THflux));

            //recolour.
            hgv2p->clearAutoscaleColour();
            if(doAutoScale){ //recolour if specified to.
                hgv1p->clearAutoscaleColour();
                hgv3p->clearAutoscaleColour();}

        }


        // rendering the hg. After each simulation step, check if enough time
        // has elapsed for it to be necessary to call v1.render().
        std::chrono::steady_clock::duration sincerender = std::chrono::steady_clock::now() - lastrender;
        if (std::chrono::duration_cast<std::chrono::milliseconds>(sincerender).count() > 17) { // 17 is about 60 Hz
            glfwPollEvents();
            v1.render();
            lastrender = std::chrono::steady_clock::now();
        }
    }

    std::cout << "Ctrl-q in graphics window to exit.\n";
    v1.keepOpen();
    return 0;
};

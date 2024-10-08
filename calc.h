#include <vector>
#include <array>
#include <sstream>
#include <iostream>
#include <morph/RD_Base.h>
#include <morph/vvec.h>
#include <functional>


/*Inherit methods and attributes from the RD_base*/
template <typename Flt>
class D_calc : public morph::RD_Base<Flt>{
public:
    alignas(alignof(std::vector<Flt>))


    //state variables.
    morph::vvec<Flt> THflux;
    morph::vvec<Flt> Fflux;
    morph::vvec<Flt> T; //temperature
    morph::vvec<Flt> show_celltype; //grid used to show rod types.

    morph::vvec<morph::vec<int,3>> coolant_positions; //RGB coordinates of coolant rods.
    morph::vvec<morph::vec<int,3>> control_positions; //RGB coordinates of coolant rods.
    morph::vvec<morph::vec<int,3>> fuel_positions; //RGB coordinates of coolant rods.
    morph::vvec<morph::vec<int,3>> source_positions; //RGB coordinates of source neutron positions.


    //coefficients & parameters.
    alignas(Flt) Flt D_Fflux = 0.1;
    alignas(Flt) Flt D_THflux = 0.1;
    alignas(Flt) Flt D_T = 0.1;

    alignas(Flt) Flt dTdt_THflux_coeff = 0.0001;
    alignas(Flt) Flt dTdt_Fflux_coeff = 0.0001;

    alignas(Flt) Flt temperature_meltdown_value = 1.0;


    alignas(Flt) Flt moderation_strength = 0.1;
    alignas(Flt) Flt absorbtion_strength = 0.1;
    alignas(Flt) Flt source_strength = 0.1;
    alignas(Flt) Flt neutrons_per_fission = 3;
    alignas(Flt) Flt fuel_strength = 1;

    alignas(bool) bool doFluxNoise = false;
    alignas(bool) bool doTemperatureNoise = false;
    alignas(bool) bool sourceNeutrons = false;
    alignas(Flt) Flt sourceStrength = 1;
    alignas(Flt) Flt noiseMaxHeight = 0;
    alignas(Flt) Flt noiseMinHeight = 0;

    //RK4 parameters.
    alignas(Flt) Flt k1 = 1.0;
    alignas(Flt) Flt k2 = 1.0;
    alignas(Flt) Flt k3 = 1.0;
    alignas(Flt) Flt k4 = 1.0;

    //get data from class D_calc is inheritting from.
    D_calc() : morph::RD_Base<Flt>() {}

        void allocate(){
        morph::RD_Base<Flt>::allocate();

        //resize the state varaiables to
        this->resize_vector_variable (this->Fflux);
        this->resize_vector_variable (this->THflux);
        this->resize_vector_variable (this->T);
        this->resize_vector_variable (this->show_celltype);
    }

    void init(){
        //zero state variables.
        this->Fflux.zero();
        this->THflux.zero();
        this->T.zero();
        this->show_celltype.zero();
        if(doFluxNoise){ //if noise is set to true, apply random values from 0 to noiseMaxHeight
            this->noiseify_vector_variable (this->Fflux, noiseMinHeight, noiseMaxHeight);
            this->noiseify_vector_variable (this->THflux, noiseMinHeight, noiseMaxHeight);
        }
        if(doTemperatureNoise){ //if noise is set to true, apply random values from noiseMinHeight to noiseMaxHeight
            this->noiseify_vector_variable (this->T, noiseMinHeight, noiseMaxHeight);
        }


        //various init function calls, see there respective functions.

        init_coolant_rods();
        init_control_rods();
        init_fuel_rods();
        init_source_positions();
        draw_celltype();
    }

    //HEX_USER_FLAG_0 IS A ***COOLANT*** channel.
    void init_coolant_rods(){
        std::list<morph::Hex>::iterator pos; //pos is a hex iterator used as a temporary variable to assign user flags to given hexes within the for loop.
        for(auto coolant_pos: this->coolant_positions){//create a control_pos vector automatically to iterate through control_postions, defined as morph::vec<int,3> above, and filled in hexvis.cpp.
            pos = this->hg->findHexAt(coolant_pos);
            if(pos != this->hg->hexen.end())
            { //if pos is the end hex, then it is probably out of range.
                pos->setUserFlags(HEX_USER_FLAG_0); //Assign the coolant channel flag to that point.
            }
        }
    }
    //HEX_USER_FLAG_1 IS A ***CONTROL*** rod.
    void init_control_rods(){
        std::list<morph::Hex>::iterator pos; //pos is a hex iterator used as a temporary variable to assign user flags to given hexes within the for loop.
        for(auto control_pos: this->control_positions){//create a control_pos vector automatically to iterate through control_postions, defined as morph::vec<int,3> above, and filled in hexvis.cpp.
            pos = this->hg->findHexAt(control_pos);
            if(pos != this->hg->hexen.end())
            { //if pos is the end hex, then it is probably out of range.
                pos->setUserFlags(HEX_USER_FLAG_1); //Assign the coolant channel flag to that point.
            }
        }
    }
    //HEX_USER_FLAG_2 IS A ***FUEL*** rod.
    void init_fuel_rods(){
        std::list<morph::Hex>::iterator pos; //pos is a hex iterator used as a temporary variable to assign user flags to given hexes within the for loop.
        for(auto fuel_pos: this->fuel_positions){ //create a fuel_pos vector automatically to iterate through fuel_postions, defined as morph::vec<int,3> above, and filled in hexvis.cpp.
            pos = this->hg->findHexAt(fuel_pos);
            if(pos != this->hg->hexen.end())
            { //if pos is the end hex, then it is probably out of range.
                pos->setUserFlags(HEX_USER_FLAG_2); //Assign the coolant channel flag to that point.
            }
        }
    }
    //HEX_USER_FLAG_3 IS A **SOURCE** position.
    void init_source_positions(){
        std::list<morph::Hex>::iterator pos; //pos is a hex iterator used as a temporary variable to assign user flags to given hexes within the for loop.
        for(auto source_pos: this->source_positions){ //create a fuel_pos vector automatically to iterate through fuel_postions, defined as morph::vec<int,3> above, and filled in hexvis.cpp.
            pos = this->hg->findHexAt(source_pos);
            if(pos != this->hg->hexen.end())
            { //if pos is the end hex, then it is probably out of range.
                pos->setUserFlags(HEX_USER_FLAG_3); //Assign the coolant channel flag to that point.
            }
        }
    }



    //compute the rate of change of fast flux
    void compute_dFfluxdt (std::vector<Flt>& Fflux_, std::vector<Flt>& dFfluxdt){
        //make a vector to then use to store laplacian
        std::vector<Flt> lapFflux(this->nhex, 0.0);
        //compute laplacian & fill the vector.
        this->compute_laplace (Fflux_, lapFflux);

        std::list<morph::Hex>::iterator hi = this->hg->hexen.begin(); //creates an iterator hi that starts on the first hex (hexen.begin)
        while(hi != this->hg->hexen.end()){ // until the end of hi, run:
            //calculating
            dFfluxdt[hi->vi] =this->D_Fflux*lapFflux[hi->vi] - this->absorbtion_strength*this->Fflux[hi->vi];

            //if the hex, hi, is a control cell.
            if(hi->getUserFlag(1)==true)
            {
                dFfluxdt[hi->vi] += -this->moderation_strength*this->Fflux[hi->vi];
            }

            //if the hex at hi is a fuel rod.
            if(hi->getUserFlag(2)==true)
            {
                dFfluxdt[hi->vi] += this->fuel_strength*(THflux[hi->vi]*this->neutrons_per_fission);
            }

            hi++; //iterate hi.
        }
    }

    //compute the rate of change of thermal flux
    void compute_dTHfluxdt (std::vector<Flt>& THflux_, std::vector<Flt>& dTHfluxdt){

        //make a vector to then use to store laplacian
        std::vector<Flt> lapTHflux(this->nhex, 0.0);
        //compute laplacian & fill the vector.
        this->compute_laplace (THflux_, lapTHflux);

        std::list<morph::Hex>::iterator hi = this->hg->hexen.begin(); //creates an iterator hi that starts on the first hex (hexen.begin)
        while(hi != this->hg->hexen.end()){ // until the end of hi, run:
            //calculating
            dTHfluxdt[hi->vi] = this->D_THflux*lapTHflux[hi->vi];

            //moderate fast flux into thermal flux if the cell is a control rod.
            if(hi->getUserFlag(1)==true)
            {
                dTHfluxdt[hi->vi] += this->moderation_strength*this->Fflux[hi->vi];
            }
            hi++; //iterate hi.

            //if the hex at hi is a fuel rod.
            if(hi->getUserFlag(2)==true)
            {
                dTHfluxdt[hi->vi] -= 0.1*THflux[hi->vi];
            }

        }
    }

    //compute rate of change of T with in input of a current state (T_) and a vector for the ROC (dTdt)
    void compute_dTdt(std::vector<Flt>& T_, std::vector<Flt>& dTdt){
        //create a vector (lapT) for & compute the scalar laplacian over the hexgrid for T_.
        std::vector<Flt> lapT(this->nhex, 0.0);
        this->compute_laplace (T_, lapT);
        std::list<morph::Hex>::iterator hi = this->hg->hexen.begin(); //creates an iterator hi that starts on the first hex, (hexen.begin)
        while(hi != this->hg->hexen.end()) // until the end of hi, run:
        {
            //calculating dT/dt
            dTdt[hi->vi] =this->dTdt_THflux_coeff*THflux[hi->vi]+this->dTdt_Fflux_coeff*Fflux[hi->vi]+this->D_T*lapT[hi->vi];
            //compute dT/dt

            //if the current hex is a coolant cell, negate a small amount of dT/dt
            //getUserFlag(0)=true if setUserFlags(HEX_USER_FLAG_0) has been run on
            //that iterator index, otherwise false. See init_coolant_rods().
            //and T at hi is greater than the subtraction I'm going to make.
            if(hi->getUserFlag(0)==true)
            {
                dTdt[hi->vi] -= 0.5*T[hi->vi];
            }
            hi++; //iterate.
        }
    }

    //adding source fast flux into a given position.
    void add_source_neutrons(){
        //pos is a hex iterator object i'm using for positioning of a hexagon.
        //Set the posittion of the neutron sources with RGB cubic coordinates. This will then be value(s) in params.json, so multiple sources can be used.

        //Add an amount of neutron flux equal to  = intial flux value * e^(-lambda * t)
        //This uses the vector index vi of pos as the location of the hexagon as a list index

        std::list<morph::Hex>::iterator hi = this->hg->hexen.begin(); //creates an iterator hi that starts on the first hex, (hexen.begin)
        while(hi != this->hg->hexen.end()) // until the end of hi, run:
        {
            if(hi->getUserFlag(3)==true)
            {
                Fflux[hi->vi] += this->source_strength*exp((-0.01)*this->stepCount);
            }
            hi++; //iterate.
        }

    }

    void draw_celltype(){ //creates the grid of values associated with the rod types.
        std::list<morph::Hex>::iterator hi = this->hg->hexen.begin(); //creates an iterator hi that starts on the first hex (hexen.begin)
        while(hi != this->hg->hexen.end()) // until the end of hi, run:
        {
            //coolant
            if(hi->getUserFlag(0)==true)
            {
                show_celltype[hi->vi] = -1;
            }
            //control
            if(hi->getUserFlag(1)==true)
            {
                show_celltype[hi->vi] = 1;
            }
            //fuel
            if(hi->getUserFlag(2)==true)
            {
                show_celltype[hi->vi] = 10;
            }


            hi++; //iterate hi.
        }
    }

    void fail_conditions(){ //check if temperature has exceeded some maimumn value
        //use a parallelized for loop to quickly go through the hexes.
        #pragma omp parallel for
        for (unsigned int h=0; h<this->nhex; ++h) { //loop from the 0th hexagon to the last hexagon.
            if(T[h]>=this->temperature_meltdown_value){ //if T at any hex exceeds a maximum, a meltdown occurs.
                //if T exceeds said value, print this message
                std::cout<<"Failure! The reactor has overheated & melt down. \n";
                std::exit(1);//then exit.
            };
        }
    }

    void step(){ //the procedure executed when the main progra steps the simulation. Calls all other relivent functions.
        if(this->sourceNeutrons){add_source_neutrons();}
        stepFflux();
        stepTHflux();
        stepT();
        fail_conditions();
        this->stepCount++; //increment stepCount.
    }

    /*
     * RK4 functions.
     *
     * They are identicle to each other. I will comment the first one 'in full'.
     */

    void stepFflux(){ //RK4 for fast flux.
        {
            // Ffluxtst: "Fflux at a test point". Ffluxtst is a temporary estimate for Fflux.
            std::vector<Flt> Ffluxtst(this->nhex, 0.0);
            std::vector<Flt> dFfluxdt(this->nhex, 0.0);
            std::vector<Flt> K1(this->nhex, 0.0);
            std::vector<Flt> K2(this->nhex, 0.0);
            std::vector<Flt> K3(this->nhex, 0.0);
            std::vector<Flt> K4(this->nhex, 0.0);
            /*
            * Stage 1
            */
            //compute the differential with the current value.
            this->compute_dFfluxdt (this->Fflux, dFfluxdt);

            #pragma omp parallel for //paralellizing the for loop.
            for (unsigned int h=0; h<this->nhex; ++h) {
                //calculating the new estimate Ffluxtst. Also calculate k1 for use later on.
                K1[h] = dFfluxdt[h] * this->dt;
                Ffluxtst[h] = this->Fflux[h] + K1[h] * 0.5 ;
            }

            /*
            * Stage 2
            */
            this->compute_dFfluxdt (Ffluxtst, dFfluxdt);
            #pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; ++h) {
                K2[h] = dFfluxdt[h] * this->dt;
                Ffluxtst[h] = this->Fflux[h] + K2[h] * 0.5;
            }

            /*
            * Stage 3
            */
            this->compute_dFfluxdt (Ffluxtst, dFfluxdt);
            #pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; ++h) {
                K3[h] = dFfluxdt[h] * this->dt;
                Ffluxtst[h] = this->Fflux[h] + K3[h];
            }

            /*
            * Stage 4
            */
            this->compute_dFfluxdt (Ffluxtst, dFfluxdt);
            #pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; ++h) {
                K4[h] = dFfluxdt[h] * this->dt;
            }

            /*
            * Final sum. This could be incorporated in the for loop for stage 4.
            */
            #pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; ++h) {
                Fflux[h] += ((K1[h] + 2.0 * (K2[h] + K3[h]) + K4[h])/(Flt)6.0);
            }
        }
    }

    void stepT(){ //RK4 for Temperature.
        {
            // Ttemp: "T at a test point". Ttemp is a temporary estimate for T.
            std::vector<Flt> Ttemp(this->nhex, 0.0);
            std::vector<Flt> dTdt(this->nhex, 0.0);
            std::vector<Flt> K1(this->nhex, 0.0);
            std::vector<Flt> K2(this->nhex, 0.0);
            std::vector<Flt> K3(this->nhex, 0.0);
            std::vector<Flt> K4(this->nhex, 0.0);
            /*
            * Stage 1
            */
            this->compute_dTdt (this->T, dTdt);
            #pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; ++h) {
                K1[h] = dTdt[h] * this->dt;
                Ttemp[h] = this->T[h] + K1[h] * 0.5 ;
            }

            /*
            * Stage 2
            */
            this->compute_dTdt (Ttemp, dTdt);
            #pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; ++h) {
                K2[h] = dTdt[h] * this->dt;
                Ttemp[h] = this->T[h] + K2[h] * 0.5;
            }

            /*
            * Stage 3
            */
            this->compute_dTdt (Ttemp, dTdt);
            #pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; ++h) {
                K3[h] = dTdt[h] * this->dt;
                Ttemp[h] = this->T[h] + K3[h];
            }

            /*
            * Stage 4
            */
            this->compute_dTdt (Ttemp, dTdt);
            #pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; ++h) {
                K4[h] = dTdt[h] * this->dt;
            }

            /*
            * Final sum. This could be incorporated in the for loop for stage 4.
            */
            #pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; ++h) {
                T[h] += ((K1[h] + 2.0 * (K2[h] + K3[h]) + K4[h])/(Flt)6.0);
            }
        }
    }

    void stepTHflux(){ //RK4 for thermal flux.
        {
            // THfluxtst: "THflux at a test point". THfluxtst is a temporary estimate for THflux.
            std::vector<Flt> THfluxtst(this->nhex, 0.0);
            std::vector<Flt> dTHfluxdt(this->nhex, 0.0);
            std::vector<Flt> K1(this->nhex, 0.0);
            std::vector<Flt> K2(this->nhex, 0.0);
            std::vector<Flt> K3(this->nhex, 0.0);
            std::vector<Flt> K4(this->nhex, 0.0);
            /*
            * Stage 1
            */
            this->compute_dTHfluxdt (this->THflux, dTHfluxdt);
            #pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; ++h) {
                K1[h] = dTHfluxdt[h] * this->dt;
                THfluxtst[h] = this->THflux[h] + K1[h] * 0.5 ;
            }

            /*
            * Stage 2
            */
            this->compute_dTHfluxdt (THfluxtst, dTHfluxdt);
            #pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; ++h) {
                K2[h] = dTHfluxdt[h] * this->dt;
                THfluxtst[h] = this->THflux[h] + K2[h] * 0.5;
            }

            /*
            * Stage 3
            */
            this->compute_dTHfluxdt (THfluxtst, dTHfluxdt);
            #pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; ++h) {
                K3[h] = dTHfluxdt[h] * this->dt;
                THfluxtst[h] = this->THflux[h] + K3[h];
            }

            /*
            * Stage 4
            */
            this->compute_dTHfluxdt (THfluxtst, dTHfluxdt);
            #pragma omp paralell for
            for (unsigned int h=0; h<this->nhex; ++h) {
                K4[h] = dTHfluxdt[h] * this->dt;
            }

            /*
            * Final sum. This could be incorporated in the for loop for stage 4.
            */
            #pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; ++h) {
                THflux[h] += ((K1[h] + 2.0 * (K2[h] + K3[h]) + K4[h])/(Flt)6.0);
            }
        }
    }

};



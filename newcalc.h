
#include <vector>
#include <array>
#include <sstream>
#include <iostream>
#include <morph/RD_Base.h>
#include <morph/vvec.h>
#include <functional>


/*Inherit methods and attributes from the RD_base*/
template <typename Flt>
class D_calc : public morph::RD_Base<Flt>
{
    public:
    alignas(alignof(std::vector<Flt>))

    morph::vvec<Flt> phi;
    morph::vvec<Flt> T; //temperature

    alignas(Flt) Flt D_phi = 0.1;
    alignas(Flt) Flt zeroPointValue = 0.0;
    alignas(bool) bool doNoise = false;
    alignas(Flt) Flt noiseHeight = 1;

    alignas(Flt) Flt k1 = 1.0;
    alignas(Flt) Flt k2 = 1.0;
    alignas(Flt) Flt k3 = 1.0;
    alignas(Flt) Flt k4 = 1.0;

    D_calc() : morph::RD_Base<Flt>() {}


    void allocate()
    {
        morph::RD_Base<Flt>::allocate();

        this->resize_vector_variable (this->phi);
        this->resize_vector_variable (this->T);
    }

    void init()
    {
        this->phi.zero();
        this->T.zero();
        if(doNoise){
        this->noiseify_vector_variable (this->phi, 0.0, noiseHeight);}
        phi[0]=zeroPointValue;
    }




    void compute_dPhidt (std::vector<Flt>& phi_, std::vector<Flt>& dPhidt)
    {
        std::vector<Flt> lapPhi(this->nhex, 0.0);
        this->compute_laplace (phi_, lapPhi);

        #pragma omp parallel for
        for (unsigned int h=0; h<this->nhex; ++h) {
            dPhidt[h] =  this->D_phi * lapPhi[h];
        }
    }



    void compute_dTdt(std::vector<Flt>& T_, std::vector<Flt>& dTdt)
    {
        std::vector<Flt> lapT(this->nhex, 0.0);
        this->compute_laplace (T_, lapT);
        #pragma omp parallel for
        for (unsigned int h=0; h<this->nhex; ++h) {
            dTdt[h] =  0.005*phi[h]+0.005*lapT[h];
        }
    }



    void step(){
        stepPhi();
        stepT();
        this->stepCount++;
    }


    void stepPhi(){
        {
            // Phitst: "Phi at a test point". Phitst is a temporary estimate for Phi.
            std::vector<Flt> Phitst(this->nhex, 0.0);
            std::vector<Flt> dPhidt(this->nhex, 0.0);
            std::vector<Flt> K1(this->nhex, 0.0);
            std::vector<Flt> K2(this->nhex, 0.0);
            std::vector<Flt> K3(this->nhex, 0.0);
            std::vector<Flt> K4(this->nhex, 0.0);
            /*
            * Stage 1
            */
            this->compute_dPhidt (this->phi, dPhidt);
            #pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; ++h) {
                K1[h] = dPhidt[h] * this->dt;
                Phitst[h] = this->phi[h] + K1[h] * 0.5 ;
            }

            /*
            * Stage 2
            */
            this->compute_dPhidt (Phitst, dPhidt);
            #pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; ++h) {
                K2[h] = dPhidt[h] * this->dt;
                Phitst[h] = this->phi[h] + K2[h] * 0.5;
            }

            /*
            * Stage 3
            */
            this->compute_dPhidt (Phitst, dPhidt);
            #pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; ++h) {
                K3[h] = dPhidt[h] * this->dt;
                Phitst[h] = this->phi[h] + K3[h];
            }

            /*
            * Stage 4
            */
            this->compute_dPhidt (Phitst, dPhidt);
            #pragma omp parallel for
                for (unsigned int h=0; h<this->nhex; ++h) {
                K4[h] = dPhidt[h] * this->dt;
            }

            /*
            * Final sum. This could be incorporated in the for loop for stage 4.
            */
            #pragma omp parallel for
            for (unsigned int h=0; h<this->nhex; ++h) {
                phi[h] += ((K1[h] + 2.0 * (K2[h] + K3[h]) + K4[h])/(Flt)6.0);
            }
        }
    }
    void stepT(){
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





};



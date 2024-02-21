
#include <vector>
#include <array>
#include <sstream>
#include <morph/RD_Base.h>
#include <morph/vvec.h>


/*Inherit methods and attributes from the RD_base*/
template <typename Flt>

class D_calc : public morph::RD_Base<Flt>{

    public:
    alignas(alignof(std::vector<Flt>))
    morph::vvec<Flt> phi;

    /*
    * The diffusion parameters.
    */
    alignas(Flt) Flt D_phi = 0.1;
    alignas(Flt) Flt zeroPointValue = 0.0;
    alignas(bool) bool doNoise = false;
    alignas(Flt) Flt noiseHeight = 1;


    /*Runge-Kutta params*/
    alignas(Flt) Flt k1 = 1.0;
    alignas(Flt) Flt k2 = 1.0;
    alignas(Flt) Flt k3 = 1.0;
    alignas(Flt) Flt k4 = 1.0;

    //Constructor
    D_calc() : morph::RD_Base<Flt>() {}

    void allocate()
    {
        // Always call allocate() from the base class first.
        morph::RD_Base<Flt>::allocate();
        // Resize and zero-initialise the various containers. Note that the size of a
        // 'vector variable' is given by the number of hexes in the hex grid which is
        // a member of this class (via its parent, RD_Base)
        this->resize_vector_variable (this->phi);
    }
    void init()
    {
        this->phi.zero();
        if(doNoise){
        this->noiseify_vector_variable (this->phi, 0.0, noiseHeight);}
        phi[0]=zeroPointValue;
    }


    /*computing dPhi/dt = D * laplacian*/
    void compute_dPhidt (std::vector<Flt>& phi_, std::vector<Flt>& dPhidt)
    {
        /*creating a vector for the laplacian over the hexgrid*/
        std::vector<Flt> lapPhi(this->nhex, 0.0);
        this->compute_laplace (phi_, lapPhi);

        #pragma omp parallel for /*for spawns a group of threads and divides loop iterations between them*/
        for (unsigned int h=0; h<this->nhex; ++h) {
            dPhidt[h] =  this->D_phi * lapPhi[h];
        }
    }



    void step()
        {
        this->stepCount++;

        // 1. 4th order Runge-Kutta computation for Phi
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
};



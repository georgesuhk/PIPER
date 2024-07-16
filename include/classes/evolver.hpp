#ifndef EVOLVER_HPP
#define EVOLVER_HPP

#include "fluxCalcs.hpp"
#include "divClean.hpp"
#include "limiter.hpp"
#include "reconstruct.hpp"

// BASE EVOLVER CLASS
class Evolver {
    public:
        Evolver(){};

        Evolver(BCFunc inputBCFunc): evolverBCFunc(inputBCFunc){};

        /* evolves the material to the next time step */
        void evolveMat(Vec2D& u, shared_ptr<SysCalcs> sysPtr, Mesh2D& mesh, double& dt, char axis);
        
        /* updates the stored flux in the evolver for the current time step */
        virtual void updateFlux(Vec2D& u, shared_ptr<SysCalcs> sysPtr, Mesh2D& mesh, double& dt, char axis);

        /* updates the fastest wave speed stored in the system for the time step */
        void updateCh(Vec2D& u, shared_ptr<SysCalcs> sysPtr, Mesh2D& mesh);

        /* returns the hyperbolic time step based on the CFL condition */
        double getTimeStep(Mesh2D& mesh);

        /* returns flux values (but mainly for debugging as flux values should not be accessed outside the evolver)*/
        Vec2D getFlux(Vec2D& u, shared_ptr<SysCalcs> sysPtr, Mesh2D& mesh, double& dt, char axis);



        // CONTROL ======

        /* returns the evolver type */
        string getEvolverType();

        /* toggles divergence cleaning */
        void setDoDC(bool inputDoDC);

        /* returns whether div clean is on */
        bool getDoDC();

    protected:

        // PARAMETERS ======

        // CFL number
        double Cnum = 0.8;

        //the hyperbolic wave speed, or the fastest speed in the system over all cells
        double ch;
        

        // CONTROL ======

        // type of evolver
        string evolverType = "base";

        // whether Divergence Cleaning is enabled
        bool doDC = false;


        // FLUXES ======

        /**
         * Note on how the flux values are stored:
         * There are nCells + 1 flux values, going from f_0.5 to f_nCells+0.5
         * These are stored in flux array index from 1 to nCells + 1
         * The 0 index of the flux array stores nothing and hence should never be accessed.
        */
        Vec2D flux;


        // BOUNDARY CONDITION ======
        BCFunc evolverBCFunc = TransBCs;
};



// FORCE (First ORder CEntered ) Scheme EVOLVER
class FORCEEvolver : public Evolver{
    public:
        FORCEEvolver(BCFunc inputBCFunc): 
        Evolver(inputBCFunc){evolverType = "FORCE"; Cnum = 0.8;};

        /* Override for how flux is updated */
        void updateFlux(Vec2D& u, shared_ptr<SysCalcs> sysPtr, Mesh2D& mesh, double& dt, char axis) override;
};



// SLIC (Slope LImited Centered) Scheme EVOLVER
class SLICEvolver : public Evolver{
    public:
        SLICEvolver(BCFunc inputBCFunc): 
        Evolver(inputBCFunc){evolverType = "SLIC"; Cnum = 0.2;};

        /* Override for how flux is updated */
        void updateFlux(Vec2D& u, shared_ptr<SysCalcs> sysPtr, Mesh2D& mesh, double& dt, char axis) override;

        void setW(double inputW);
        double getW();

        void setLimiter(Limiter inputLimFunc);
        Limiter getLimiter();

    protected:
        double w = 1;
        Limiter limFunc = vanLeer;
};












#endif
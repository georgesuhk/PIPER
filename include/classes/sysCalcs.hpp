#ifndef SYSCALC_HPP
#define SYSCALC_HPP

#include "settings.hpp"
#include "EoS.hpp"


//base class for system calculations
class SysCalcs {
    public:
        SysCalcs(shared_ptr<EoS> inputEoSPtr): EoSPtr(inputEoSPtr){};


        // STATE VARIABLE EVALUATION & CONVERSIONS ======
        // Operated on primitive variables -------
        
        /* returns internal energy */
        double get_e_Prim(CellVec& u);

        /* returns the temperature */
        double get_T_Prim(CellVec& uPrim);

        /* returns total non-magnetic energy (kinetic + internal) */
        double get_E_Prim(CellVec& uPrim);

        /* returns total energy (kinetic + internal + magnetic) */
        double get_U(CellVec& uPrim);

        /* returns bulk flow kinetic energy*/
        double get_KE_Prim(CellVec& uPrim);

        /* converts CellVec from primitive form to conservative form */
        virtual CellVec primToConserv(CellVec& uPrim) = 0;


        // Operated on conservative variables -------
        /* returns the temperature*/
        double get_T(CellVec& u);

        /* returns total non-magnetic energy (kinetic + internal) */
        double get_E(CellVec& u);

        /* returns internal energy */
        double get_e(CellVec& u);

        /* returns pressure through interpolation */ 
        double get_p(CellVec& u);  

        /* returns pressure via cached values */
        double get_p(int i, int j);

        /* returns bulk flow kinetic energy*/
        double get_KE(CellVec& u);

        /* returns the resistivity */
        double get_Resis(CellVec& u, int i, int j, bool interp = true);

        /* converts CellVec from conservative form to primitive form */
        virtual CellVec conservToPrim(CellVec& u) = 0;


        // Can be operated on both ------
        double get_MagE(CellVec& u);



        // FLUX FUNCTIONS =====
        virtual CellVec f(CellVec& u) = 0;

        virtual CellVec g(CellVec& u) = 0;


        
        // WAVE SPEED CALCULATIONS =====
        /* returns sound speed */
        double get_Cs(CellVec& u);

        /* returns Alfven wave speed */
        double get_Ca(CellVec& u, char axis);

        /* returns fast magnetosonic speed */
        double get_Cf(CellVec& u, char axis);

        /* returns the speed of the fastest wave in the system (hyperbolic wave speed) */
        double getFastestWaveSpeed(CellVec& u);



        // OTHER EVALUATIONS =====
        Scalar2D getDivBField(Vec2D& u, Mesh2D& mesh);



        // OTHER ======
        string getSysName();
        shared_ptr<EoS> getEoSPtr();

    protected:
        shared_ptr<EoS> EoSPtr;
        string sysName = "base";

};

// CALCULATIONS FOR A FULLY IONISED PLASMA
class FIPCalcs : public SysCalcs{
    public:
        FIPCalcs(shared_ptr<EoS> inputEoSPtr):SysCalcs(inputEoSPtr){sysName = "FIP";};


        // STATE VARIABLE EVALUATION & CONVERSIONS ======
        // operated on primitive variables ------
        virtual CellVec primToConserv(CellVec& uPrim) override;

        // Operated on conservative variables -------
        virtual CellVec conservToPrim(CellVec& u) override;


        // FLUX FUNCTIONS ======
        virtual CellVec f(CellVec& u) override;
        virtual CellVec g(CellVec& u) override;

    protected:
        int cellVecLen = 9; // number of variables stored for a cell, local to this sysCalc
      
};

#endif
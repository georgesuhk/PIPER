#ifndef SYSCALC_HPP
#define SYSCALC_HPP

#include "settings.hpp"
#include "EoS.hpp"
#include "calcs.hpp"
#include "sysEigens.hpp"
// #include <Eigen/Dense>

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

        /* returns mass fraction of ions */
        double interp_mass_frac_i_Prim(CellVec& uPrim);


        /* converts CellVec from primitive form to conservative form */
        virtual CellVec primToConserv(CellVec& uPrim) = 0;


        // Operated on conservative variables -------

        /* returns the temperature*/
        double get_T(int i, int j);

        double interp_T(CellVec& u);
        
        double interp_T(double& rho, double& p);

        /* returns total non-magnetic energy (kinetic + internal) */
        double get_E(CellVec& u);

        /* returns internal energy */
        double get_e(CellVec& u);

        /* returns the adibatic index / ratio between specific heats */
        double interp_gamma(CellVec& u);

        /* returns thermal conductivity */
        double interp_therm_con(CellVec& u);

        /* returns pressure through interpolation */ 
        double interp_p(CellVec& u);  

        /* returns pressure via cached values */
        double get_p(int i, int j);

        /* returns bulk flow kinetic energy*/
        double get_KE(CellVec& u);

        /* calculates KE based on given velocities */
        double get_KE(double& rho, double& vx, double& vy, double& vz);

        /* calculate KE based on velocity squared */
        double get_KE(double& rho, double v_squared);

        /* returns velocity squared */
        double get_v_squared(CellVec& u);

        /* returns the resistivity */
        double get_Resis(CellVec& u, int i, int j, bool interp = true);

        double interp_Resis(double& rho, double& p);

        double interp_Resis(CellVec& u);

        // Get masses

        /* returns neutral mass */
        double get_m_n();

        /* returns ion mass */
        double get_m_i();

        // Species densities ------

        /* electrons */
        double interp_mass_frac_e(CellVec& u);

        /* ions */        
        double interp_mass_frac_i(CellVec& u);

        double interp_mass_frac_i(double& rho, double& p);

        double get_mass_frac_i(int i, int j);

        double interp_e_n(CellVec& u);

        double interp_e_n(double& rho, double& p);

        /* converts CellVec from conservative form to primitive form */
        virtual CellVec conservToPrim(CellVec& u) = 0;


        // Can be operated on both ------
        double get_MagE(CellVec& u);



        // FLUX FUNCTIONS =====
        virtual CellVec f(CellVec& u, int i, int j, bool interp = false) = 0;

        virtual CellVec g(CellVec& u, int i, int j, bool interp = false) = 0;


        
        // WAVE SPEED CALCULATIONS =====
        /* returns sound speed */
        double get_Cs(CellVec& u);

        /* returns Alfven wave speed */
        double get_Ca(CellVec& u, char axis);

        /* returns fast magnetosonic speed */
        double get_Cf(CellVec& u, char axis);

        /* returns the speed of the fastest wave in the system (hyperbolic wave speed) */
        virtual double getFastestWaveSpeed(CellVec& u) = 0;



        // OTHER EVALUATIONS =====
        Scalar2D getDivBField(Vec2D& u, Mesh2D& mesh);

        /* calculates and sets w within u based on the no inertial approximation */
        virtual void set_w_no_inert(Vec2D& u, Mesh2D& mesh) = 0;



        // OTHER ======
        string getSysName();
        shared_ptr<EoS> getEoSPtr();

    protected:
        shared_ptr<EoS> EoSPtr;
        string sysName = "base";
        double m_i = protonMass;
        double m_n = protonMass; 
        Eigen::Matrix<double, 7, 7> eigen_matrix;
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
        virtual CellVec f(CellVec& u, int i, int j, bool interp = false) override;
        virtual CellVec g(CellVec& u, int i, int j, bool interp = false) override;


        // WAVE CALCULATIONS ======
        virtual double getFastestWaveSpeed(CellVec& u) override;


    protected:
        int cellVecLen = 9; // number of variables stored for a cell, local to this sysCalc
      
};


// CALCULATIONS FOR A PARTIALLY IONISED PLASMA

/**
 * PIP_0:
 * Frame for PIP systems
 * Non-conductive
 * Ignores inertia in dw/dt, i.e. uses w is not evolved explicitly
 */
class PIP0_Calcs : public SysCalcs{
    public:
        PIP0_Calcs(shared_ptr<EoS> inputEoSPtr):SysCalcs(inputEoSPtr){sysName = "PIP0";};


        // STATE VARIABLE EVALUATION & CONVERSIONS ======
        // operated on primitive variables ------
        virtual CellVec primToConserv(CellVec& uPrim) override;

        // Operated on conservative variables -------
        virtual CellVec conservToPrim(CellVec& u) override;
        double get_vix(CellVec& u);

        /* returns the kinectic + internal energy for neutrals only */
        double get_total_neutral_E(CellVec& u);


        // FLUX FUNCTIONS ======
        virtual CellVec f(CellVec& u, int i, int j, bool interp = false) override;
        virtual CellVec g(CellVec& u, int i, int j, bool interp = false) override;


        // WAVE CALCULATIONS ======
        virtual double getFastestWaveSpeed(CellVec& u) override;


        // OTHER ======

        /* calculates and sets w within u based on the no inertial approximation */
        virtual void set_w_no_inert(Vec2D& u, Mesh2D& mesh) override;

    protected:
        int cellVecLen = 12; // number of variables stored for a cell, local to this sysCalc   
};



// // TWO FLUID PARTIALLY IONISED PLASMA

// /**
//  * PIP_0:
//  * Frame for PIP systems
//  * Non-conductive
//  * Ignores inertia in dw/dt, i.e. uses w is not evolved explicitly
//  */
// class PIP2F_Calcs : public SysCalcs{
//     public:
//         PIP2F_Calcs(shared_ptr<EoS> inputEoSPtr):SysCalcs(inputEoSPtr){sysName = "PIP2F";};


//         // STATE VARIABLE EVALUATION & CONVERSIONS ======
//         // operated on primitive variables ------
//         virtual CellVec primToConserv(CellVec& uPrim) override;
//         double get_KE_Prim_neutral(CellVec& uPrim);
//         double get_E_Prim_neutral(CellVec& uPrim);

//         // Operated on conservative variables -------
//         virtual CellVec conservToPrim(CellVec& u) override;
//         array<double,2> get_p_2F(CellVec& u);

//         // FLUX FUNCTIONS ======
//         virtual CellVec f(CellVec& u, int i, int j, bool interp = false) override;
//         virtual CellVec g(CellVec& u, int i, int j, bool interp = false) override;

//         // WAVE CALCULATIONS ======
//         virtual double getFastestWaveSpeed(CellVec& u) override;


//     protected:
//         int cellVecLen = 14; // number of variables stored for a cell, local to this sysCalc   
// };

#endif
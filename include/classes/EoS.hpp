#ifndef EOS_HPP
#define EOS_HPP

#include "settings.hpp"
#include "utils.hpp"


// Base Class
class EoS {
    public:
        /**
         * Default (base) constructor
         * By calling EoS with a default constructor, a EoS object of "base" type will be created that should not be
         * used in simulation, and acts to signal that the EoS has not yet been proper initialised.
        */
        EoS():EoSType("base"){};

        string getEoSType();

        /* returns specific internal energy */
        virtual double get_e(double& rho, double& p) = 0;

        /* returns temperature */
        virtual double get_T(double& rho, double& p) = 0;

        /* returns pressure */
        virtual double get_p(double& rho, double& e, bool interp) = 0;

        /* returns sound speed */
        virtual double get_Cs(double& rho, double& p) = 0;

        /* returns resistivity*/
        virtual double get_Resis(int& i, int& j) = 0;
        virtual double interp_Resis(double& rho, double& p);

    
    protected:
        string EoSType;
    };

// Ideal EoS
class IdealEoS : public EoS {
    public:
        IdealEoS(double inputGamma, double input_m): 
        gamma(inputGamma), m(input_m){EoSType = "Ideal";};

        // GET FUNCTIONS ======
        /* returns specific internal energy */
        virtual double get_e(double& rho, double& p) override;

        /* returns temperature */
        virtual double get_T(double& rho, double& p) override;

        /* returns pressure */
        virtual double get_p(double& rho, double& e, bool interp) override;

        /* returns sound speed */
        virtual double get_Cs(double& rho, double& p) override;

        /* returns resistivity*/
        virtual double get_Resis(int& i, int& j) override;

        double get_gamma();

        // SETTING VARIABLES ======
        void set_gamma(double inputGamma);
        void set_constResis(double inputResis);
    
    protected:
        // adibatic index
        double gamma;
        // (average) particle mass
        double m;
        // resistivity
        double constResis = 0;
    };

class TabEoS : public EoS{
    public: 
        TabEoS():EoS(){EoSType = "Tab";};

        // GET FUNCTIONS ======
        /* returns the length of the pressure idx */
        int get_pLen();

        /* returns the length of the density idx */
        int get_rhoLen();

        /* returns specific internal energy */
        virtual double get_e(double& rho, double& p) override;

        /* returns temperature */
        virtual double get_T(double& rho, double& p) override;

        /* returns pressure */
        virtual double get_p(double& rho, double& e, bool interp) override;

        /* returns sound speed */
        virtual double get_Cs(double& rho, double& p) override;

        /* returns resistivity*/
        virtual double get_Resis(int& i, int& j) override; // note properly implemented yet

        /* interpolates resistivity based on density and pressure */
        virtual double interp_Resis(double& rho, double& p) override;

        // CACHING ======
        // interpVars function

        // GENERATION OF TABEOS ======
        void genFromData(Mesh2D mesh, vector<string> varList, string dataFolder, char delimiter);


    protected:
        // array-like data
        Scalar1D densities;
        Scalar1D pressures;  

        // tabular data
        Scalar2D e_Table;
        Scalar2D T_Table;
        Scalar2D Cs_Table;
        Scalar2D elecConduct_Table;
        Scalar2D thermConduct_Table;

        // Caching
        Scalar2D resis_Cache;

        // Control variables
        bool cached = false;
        array<int,2> activeRhoIndices;
        array<int,2> activePIndices;

        


};

#endif
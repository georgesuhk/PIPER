#ifndef EOS_HPP
#define EOS_HPP

#include "settings.hpp"
#include "utils.hpp"
#include "parse.hpp"


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
        virtual double get_p(int& i, int& j) = 0;

        virtual double interp_p(double& rho, double& e) = 0;
 
        /* returns sound speed */
        virtual double get_Cs(double& rho, double& p) = 0;

        /* returns resistivity */
        virtual double get_Resis(int& i, int& j) = 0;

        /* interp resistivity */
        virtual double interp_Resis(double& rho, double& p);

        // CACHING

        /* caches all cached variables */
        void cacheAll(Vec2D& u, Mesh2D& mesh);

        /* caches pressure */
        void cache_p(Vec2D& u, Mesh2D& mesh);

        /* caches resistivity */
        void cache_resis(Vec2D& u, Mesh2D& mesh);

    
    protected:
        string EoSType;

        // Caching ------
        Scalar2D e_cache;
        Scalar2D p_cache;
        Scalar2D resis_cache;

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
        virtual double get_p(int& i, int& j) override;

        /* interps pressure */
        virtual double interp_p(double& p, double& e) override;

        /* returns sound speed */
        virtual double get_Cs(double& rho, double& p) override;

        /* returns resistivity */
        virtual double get_Resis(int& i, int& j) override;

        /* interps resistivity */
        virtual double interp_Resis(double& rho, double& p) override;

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

        /* returns the entire pressures array */
        Scalar1D& get_pressures();

        /* returns the entire densities array */
        Scalar1D& get_densities();

        /* returns the entire temperatures table */
        Scalar2D& get_T_Table();

        /* returns specific internal energy */
        virtual double get_e(double& rho, double& p) override;

        /* returns temperature */
        virtual double get_T(double& rho, double& p) override;

        /* returns pressure */
        virtual double get_p(int& i, int& j) override;

        /* interps p */
        virtual double interp_p(double& rho, double& e) override;

        /* returns sound speed */
        virtual double get_Cs(double& rho, double& p) override;

        /* returns resistivity*/
        virtual double get_Resis(int& i, int& j) override;

        /* interpolates resistivity based on density and pressure */
        virtual double interp_Resis(double& rho, double& p) override;

        // CACHING ======
        // interpVars function

        // GENERATION OF TABEOS ======
        void genFromData(Mesh2D mesh, vector<string> varList, string dataFolder, char delimiter);


    protected:
        // array-like data ------
        Scalar1D densities;
        Scalar1D pressures;  

        // tabular data ------

        /* coefficients */
        Scalar2D e_Table;
        Scalar2D T_Table;
        Scalar2D Cs_Table;
        Scalar2D resis_Table;
        Scalar2D thermConduct_Table;

        /* species densities */
        Scalar2D n_e_Table;
        Scalar2D n_n_Table;
        Scalar2D n_i_Table;


        // Control variables ------
        array<int,2> activeRhoIndices;
        array<int,2> activePIndices;
};

#endif
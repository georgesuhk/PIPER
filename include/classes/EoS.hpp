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

        /* returns ratio of specific heats */
        virtual double interp_gamma(double& rho, double& p) = 0;

        virtual double get_gamma(int& i, int& j) = 0;

        /* returns temperature */
        virtual double interp_T(double& rho, double& p) = 0;

        virtual double get_T(int& i, int& j) = 0;

        /* returns pressure */
        virtual double interp_p(double& rho, double& e) = 0;

        virtual double get_p(int& i, int& j) = 0;

        /* returns sound speed */
        virtual double get_Cs(double& rho, double& p) = 0;

        /* returns mass fraction of electrons */
        virtual double interp_mass_frac_e(double& rho, double& p) = 0;

        /* returns mass fraction of ions */
        virtual double interp_mass_frac_i(double& rho, double& p) = 0;

        virtual double get_mass_frac_i(int& i, int& j) = 0;

        /* returns resistivity */
        virtual double interp_Resis(double& rho, double& p);

        virtual double get_Resis(int& i, int& j) = 0;

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
        Scalar2D gamma_cache;
        Scalar2D e_cache;
        Scalar2D p_cache;
        Scalar2D resis_cache;
        Scalar2D mfi_cache;
        Scalar2D T_cache;

    };

// Ideal EoS
class IdealEoS : public EoS {
    public:
        IdealEoS(double inputGamma, double input_m): 
        gamma(inputGamma), m(input_m){EoSType = "Ideal";};

        // GET FUNCTIONS ======

        /* returns specific internal energy */
        virtual double get_e(double& rho, double& p) override;

        /* returns adiabatic index / ratio of specific heats through interpolation */
        double interp_gamma(double& rho, double& p) override;

        /* returns adiabatic index from cached values */
        double get_gamma(int& i, int& j) override;

        /* returns temperature */
        virtual double interp_T(double& rho, double& p) override;

        virtual double get_T(int& i, int& j) override;

        /* returns pressure */
        virtual double interp_p(double& p, double& e) override;

        virtual double get_p(int& i, int& j) override;

        /* returns sound speed */
        virtual double get_Cs(double& rho, double& p) override;

        /* returns mfe */
        virtual double interp_mass_frac_e(double& rho, double& p) override;

        /* returns mfi */
        virtual double interp_mass_frac_i(double& rho, double& p) override;

        virtual double get_mass_frac_i(int& i, int& j) override;

        /* returns resistivity */
        virtual double get_Resis(int& i, int& j) override;

        /* interps resistivity */
        virtual double interp_Resis(double& rho, double& p) override;

        // SETTING VARIABLES ======
        void set_gamma(double inputGamma);
        void set_constResis(double inputResis);
        void set_mass_frac_n(double input_mass_frac);
        void set_mass_frac_i(double input_mass_frac);
    
    protected:
        // adibatic index
        double gamma;
        // (average) particle mass
        double m;
        // resistivity
        double constResis = 0;
        // constant mass fracs
        double mass_frac_e = 0;
        double mass_frac_i = 1.0;
        double mass_frac_n = 0.0;
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

        /* returns adiabatic index / ratio of specific heats */
        double interp_gamma(double& rho, double& p) override;

        double get_gamma(int& i, int& j) override;

        /* returns temperature */
        virtual double interp_T(double& rho, double& p) override;

        virtual double get_T(int& i, int& j) override;

        /* returns pressure */
        virtual double get_p(int& i, int& j) override;

        /* interps p */
        virtual double interp_p(double& rho, double& e) override;

        /* returns sound speed */
        virtual double get_Cs(double& rho, double& p) override;

        /* returns mfe */
        virtual double interp_mass_frac_e(double& rho, double& p) override;

        /* returns mfi */
        virtual double interp_mass_frac_i(double& rho, double& p) override;

        virtual double get_mass_frac_i(int& i, int& j) override;

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
        Scalar2D gamma_Table;
        Scalar2D T_Table;
        Scalar2D Cs_Table;
        Scalar2D resis_Table;
        Scalar2D thermConduct_Table;

        /* species specific data */
        Scalar2D mass_frac_e_Table;
        Scalar2D mass_frac_n_Table;
        Scalar2D mass_frac_i_Table;


        // Control variables ------
        array<int,2> activeRhoIndices;
        array<int,2> activePIndices;
};

#endif
#ifndef CALCS_HPP
#define CALCS_HPP

#include "constants.hpp"
// GENERAL CALCS ======

/* returns temperature according to ideal gas equation */
double getIdealGasTemp(double rho, double p, double particleMass);

/* returns density according to ideal gas equation */
double getIdealGasDensity(double p, double temp, double particleMass);

/* returns pressure according to the ideal gas equation */
double getIdealGasPressure(double rho, double temp, double particleMass);

/* returns the standard reduced mass */
double get_reduced_mass(double m1, double m2);

/* converts temperature from Kelvin to eV */
double K_to_eV(double T_K);

/* returns the Coulomb logarithm */
double getCoulombLog(double& n_e, double T);

/* returns the number density of ions */
double get_n_i(double rho, double mass_frac_i, double m_i);

/* returns the number density of neutrals */
double get_n_n(double rho, double mass_frac_n, double m_n);

/* returns the number density of electrons */
double get_n_e(double rho, double mass_frac_e);



// COLLISION RATES (nu) ======

/* returns nu_in, the collision rate between ions and neutrals */
double get_coll_freq_in(double& n_n, double& m_i, double& m_n, double& T);

/* returns nu_en, the collision rate between e- and neutrals */
double get_coll_freq_en(double n_n, double& m_i, double& m_n, double& T);

/* returns nu_ei, the collision rate between e- and ions */
double get_coll_freq_ei(double& n_i, double& coulombLog, double& T);

// COLLISION COEFFICIENTS (alpha) ======

/* returns the i-n alpha collision coefficients found in Braginskii 
   * defined such that R_in = - alpha_in * (vi - vn)
*/
double get_coll_coeff_in(double& n_i, double& n_n, double& m_i, double m_n, double& T);

/* returns the e-n alpha collision coefficients found in Braginskii 
   * defined such that R_en = - alpha_en * (ve - vn)
*/
double get_coll_coeff_en(double& n_i, double& n_n, double& m_i, double& m_n, double& T);

/* returns the e-i alpha collision coefficient */
double get_coll_coeff_ei(double n_e, double& n_i, double& coulombLog, double& T);

// non scaled versions

double get_coll_freq_in_noScale(double& n_n, double& m_i, double& m_n, double& T);

double get_coll_freq_en_noScale(double n_n, double& m_i, double& m_n, double& T);

double get_coll_freq_ei_noScale(double& n_i, double& coulombLog, double& T);

double get_coll_coeff_in_noScale(double& n_i, double& n_n, double& m_i, double m_n, double& T);

double get_coll_coeff_en_noScale(double& n_i, double& n_n, double& m_i, double& m_n, double& T);

double get_coll_coeff_ei_noScale(double n_e, double& n_i, double& coulombLog, double& T);




#endif
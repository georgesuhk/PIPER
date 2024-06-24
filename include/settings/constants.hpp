#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include "modules.hpp"

// constants
const double pAtmos = 101325.0;
const double myPI = 3.14159265358979323846;
const double kB = 1.38e-23;
const double planckConst = 6.626e-34;
const double dalton = 1.661e-27;
const double uniGasConst = 8.314;
const double numAvogadro = 6.022e23;
const double elementaryCharge = 1.602e-19;
const double electronMass = 9.109e-31;
const double protonMass = 1.6726e-27;
const double vacPermeab = 1.2566e-6;

// material properties
const double airParticleMass = 4.81e-26;
const double hydrogenMass = 1.00784 * dalton;
const double H_ionize_energy = 13.6 * elementaryCharge;

// coefficients
const double sigma_coll_in = 5e-19;
const double sigma_coll_en = 1e-19;

/* coefficients from Braginskii 1965 */
const double brag_a0 = 0.5063;
// const double brag_a0 = 1;


// scaling coefficients

/* scaling factor for Amps and seconds (due to scaling for p and B) -> to put system in atms for p and B/sqrt(u0) for B */
const double resisScaling = 1 / (vacPermeab * sqrt(pAtmos));




#endif
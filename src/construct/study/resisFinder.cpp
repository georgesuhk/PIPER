#include "mutation++.h"
#include "EoSUtils.hpp"
#include "parse.hpp"

using namespace Mutation;
using namespace Mutation::Thermodynamics;
using namespace Mutation::Transport; 
 
const int cellVarsNums = 12;
const int omp_threads = 4;

// SETTINGS ------

double p_SF = 1;
double rho_SF = 1;

char delimiter = ' ';
double particleMass = protonMass;
string mixName = "HFusion";

double p = 101325 / 1e4;
// double p = 10;
double rho = 1.67e-7;
double resisScale = 1.0/(sqrt(pAtmos)*vacPermeab);

int main(void){


    Mixture mix(mixName);

    double temp = BisectSolverEoS_T(rho, p, particleMass, mix, 0.01, 1e-9, 10000, 0.5);

    mix.equilibrate(temp, p);

    double number_density = mix.numberDensity();

    double mass_frac_e = mix.Y()[0];
    double mass_frac_n = mix.Y()[1];
    double mass_frac_i = mix.Y()[2];

    double density = mix.density();
    double n_e = number_density * mix.X()[0];
    double n_n = number_density * mix.X()[1];
    double n_i = number_density * mix.X()[2];

    double cLog = getCoulombLog(n_e, temp);

    double coll_freq_in = get_coll_freq_in_noScale(n_n, particleMass, particleMass, temp);
    double coll_freq_en = get_coll_freq_en_noScale(n_n, particleMass, particleMass, temp);
    double coll_freq_ei = get_coll_freq_ei_noScale(n_i, cLog, temp);

    double alpha_in = get_coll_coeff_in_noScale(n_i, n_n, particleMass, particleMass, temp);
    double alpha_en = get_coll_coeff_en_noScale(n_i, n_n, particleMass, particleMass, temp);
    double alpha_ei = get_coll_coeff_ei_noScale(n_e, n_i, cLog, temp);

    double gamma = mix.mixtureEquilibriumGamma();
    double resis_mpp = 1/mix.electricConductivity();
    double resis_Brag = (alpha_en + alpha_ei) / (pow(n_e, 2.0) * pow(elementaryCharge, 2.0));
    double e = mix.mixtureEnergyMass();
    vector<double> eArray(3);
    mix.getEnergiesMass(eArray.data());


    cout << "T: " << temp << endl;
    cout << "P: " << p << endl;
    cout << "e: " << e << endl;
    cout << "cv * T: " << mix.mixtureFrozenCvMole() * temp << endl;
    cout << "cp * T: " << mix.mixtureFrozenCpMole() * temp << endl;
    // cout << "cv_eq * T: " << mix.mixtureEquilibriumCvMass() * temp << endl;
    cout << "cp_eq * T: " << mix.mixtureEquilibriumCpMole() * temp << endl;



    cout << "rho * e: " << rho * e << endl;
    cout << "eArray sum: " << (mass_frac_e * eArray[0] + mass_frac_n * eArray[1] + mass_frac_i * eArray[2]) * rho << endl;

    cout << "rho_n e_n: " << mass_frac_n * rho * eArray[1] << endl;
    cout << "mfn * rho * e: " << mass_frac_n * rho * e << endl;

    // cout << "e: " << e << endl;
    double effective_mass_fac = (mass_frac_e * electronMass + mass_frac_i * protonMass + mass_frac_n * protonMass) / protonMass;
    cout << "(gamma)/(gamma-1) p: " << (gamma / (gamma - 1)) * p * mix.mixtureMwMole(mix.X()) << endl;

    cout << "e_n = " << mix.X()[1] * e << endl;
    cout << "eArray n: " << mass_frac_n * rho * eArray[1] << endl;

    cout << "n kB T: " << number_density * kB * temp << endl;
    cout << "approx p: " << (n_n + 2 * n_i) * kB * temp << endl;
    cout << "approx e: " << ((n_n + 2 * n_i) * kB * temp) / (gamma - 1) << endl;


    cout << "rho: " << rho << endl;
    cout << "n: " << number_density << endl;
    cout << "mfi: " << mass_frac_i << endl;
    cout << "mfn: " << mass_frac_n << endl;
    cout << "gamma: " << gamma << endl;

    cout << "\n" << endl;
    cout << "resis_mpp: " << resis_mpp << "; scaled: " << resis_mpp * resisScale << endl;
    cout << "resis_Brag: " << resis_Brag << "; scaled: " << resis_Brag * resisScale << endl;

    vector<double> enthalpies(3);
    mix.getEnthalpiesMass(enthalpies.data());

    cout << "h_n: " << enthalpies[2] - n_i * kB * temp / (mass_frac_i * rho) << endl;
    cout << "e_n: " << eArray[2] << endl;

    /* testing set state */


    vector<double> eTotalArray = {density * e};

    vector<double> densityArray = {density*(mass_frac_i/1836.15), density*mass_frac_n, density*mass_frac_i};
    vector<double> densityArray2 = {density/1000};

    vector<double> tempArray = {temp};

    
    
    mix.setState(densityArray.data(), tempArray.data(), 1);
    cout << "pressure from set state: " << mix.P() << endl;

    mix.setState(densityArray2.data(), eTotalArray.data(), 0);
    cout << "pressure from set state: " << mix.P() << endl;
    
    return 0; 
}
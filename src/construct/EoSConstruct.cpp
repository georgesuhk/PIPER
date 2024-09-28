#include "mutation++.h"
#include "simulation.hpp"
#include "EoSUtils.hpp"
#include "parse.hpp"

using namespace Mutation;
using namespace Mutation::Thermodynamics;
using namespace Mutation::Transport; 
 
const int cellVarsNums = 9;
const int omp_threads = 1;
double rho_SF = 1;
double p_SF = 1;

// SETTINGS ------

char delimiter = ' ';
string folder = "./EoSData/HFusion4_S4/";
double particleMass = protonMass;
string mixName = "HFusion";

// SCALINGS ------

double pScale = 1.0/pAtmos;
double eScale = 1.0/pAtmos;
double CsScale = 1.0/sqrt(pAtmos);
double resisScale = 1.0/(sqrt(pAtmos)*vacPermeab);
double thermConScale = 1.0/(pow(pAtmos,1.5));

// ROUNDING ------
double roundFactor = 1.0e8;

int main(void){

    // LOADING IN P AND RHO RANGES ======

    vector<double> pressures = singleColParse(folder + "pRange.csv");
    vector<double> densities = singleColParse(folder + "rhoRange.csv");

    int pLen = pressures.size();
    int rhoLen = densities.size();

    cout << "Constructing EoS, folder: " << folder << endl;
    cout << "pLen: " << pLen << ", rhoLen: " << rhoLen << endl;

    // CREATING STORAGE FOR ALL VARIABLES ======
    Mixture mix("HFusion");

    Scalar2D pIndex = makeScalar2D(pLen, rhoLen);
    Scalar2D rhoIndex = makeScalar2D(pLen, rhoLen);

    Scalar2D T_Table = makeScalar2D(pLen, rhoLen);
    Scalar2D e_Table = makeScalar2D(pLen, rhoLen);
    Scalar2D e_n_Table = makeScalar2D(pLen, rhoLen); // specific internal energy for the neutrals
    Scalar2D gamma_Table = makeScalar2D(pLen, rhoLen);
    Scalar2D Cs_Table = makeScalar2D(pLen, rhoLen);
    Scalar2D resis_Table = makeScalar2D(pLen, rhoLen);
    Scalar2D thermCon_Table = makeScalar2D(pLen, rhoLen);  

    Scalar2D mass_frac_e_Table = makeScalar2D(pLen, rhoLen);  
    Scalar2D mass_frac_n_Table = makeScalar2D(pLen, rhoLen);   
    Scalar2D mass_frac_i_Table = makeScalar2D(pLen, rhoLen);   

    // CONSTRUCTING ======
    for (int pIdx = 0; pIdx < pressures.size(); pIdx++){
        cout << "\rProgress: " << pIdx << "/" << pressures.size() << "     " <<std::flush;
        for (int rhoIdx = 0; rhoIdx < densities.size(); rhoIdx++){

            double pressure = pressures[pIdx];
            double rho = densities[rhoIdx];  
            double temp = BisectSolverEoS_T(rho, pressure, particleMass, mix, 1, 0.001, 200, 0.5);
            
            pIndex[pIdx][rhoIdx] = pressure * pScale;
            rhoIndex[pIdx][rhoIdx] = rho;
            T_Table[pIdx][rhoIdx] = temp;

            // EQUILIBRATE SYSTEM AND GET VALUES FROM Mutation++ ------
            mix.equilibrate(temp, pressure);

            /* Coefficients -> needs to be scaled */ 

            e_Table[pIdx][rhoIdx] = mix.mixtureEnergyMass() * eScale;
            vector<double> eArray(3);
            mix.getEnergiesMass(eArray.data());
            e_n_Table[pIdx][rhoIdx] = eArray[1] * eScale;
            gamma_Table[pIdx][rhoIdx] = mix.mixtureEquilibriumGamma();
            Cs_Table[pIdx][rhoIdx] = mix.equilibriumSoundSpeed() * CsScale;
            
            double n_e = mix.numberDensity() * mix.X()[0];
            double n_n = mix.numberDensity() * mix.X()[1];
            double n_i = mix.numberDensity() * mix.X()[2];
            double cLog = getCoulombLog(n_e, temp);

            double alpha_in = get_coll_coeff_in_noScale(n_i, n_n, particleMass, particleMass, temp);
            double alpha_en = get_coll_coeff_en_noScale(n_i, n_n, particleMass, particleMass, temp);
            double alpha_ei = get_coll_coeff_ei_noScale(n_e, n_i, cLog, temp);
            double coll_freq_in = get_coll_freq_in_noScale(n_n, particleMass, particleMass, temp);
            double coll_freq_ei = get_coll_freq_ei_noScale(n_i, cLog, temp);

            /* resistivity */
            // double resis = 1/mix.electricConductivity();
            double resis = (alpha_en + alpha_ei) / (pow(n_e, 2.0) * pow(elementaryCharge, 2.0));
            resis_Table[pIdx][rhoIdx] = resis * resisScale;

            /* thermal conductivity */
            // double thermCon = mix.equilibriumThermalConductivity();
            // cout << "temp: " << temp << endl;

            double thermCon_n = 2.4 * (n_n * kB * temp)/(particleMass * coll_freq_in);
            // cout << "n: " << thermCon_n << endl;
        
            double thermCon_i = 3.906 * (n_i * kB * temp) / (particleMass * coll_freq_ei);
            // cout << "i: " << thermCon_i << endl;

            double thermCon_e = 3.16 * (n_e * kB * temp) / (electronMass * coll_freq_ei);
            // cout << "e: " << thermCon_e << endl;
            
            double thermCon = thermCon_n + thermCon_i + thermCon_e;
            // cout << "therm Con Scale: " << thermConScale << endl;

            thermCon_Table[pIdx][rhoIdx] = thermCon * thermConScale;


            /* Number densities */
            mass_frac_e_Table[pIdx][rhoIdx] = round(mix.Y()[0] * roundFactor) / roundFactor;
            mass_frac_n_Table[pIdx][rhoIdx] = round(mix.Y()[1] * roundFactor) / roundFactor;
            mass_frac_i_Table[pIdx][rhoIdx] = round(mix.Y()[2] * roundFactor) / roundFactor;

        }
    }

    cout << "Exporting EoS" << endl; 

    easyExport(pIndex, "pressure", folder);
    easyExport(rhoIndex, "densities", folder);

    easyExport(T_Table, "T", folder);
    easyExport(e_Table, "e", folder);
    easyExport(e_n_Table, "e_n", folder);
    easyExport(gamma_Table, "gamma", folder);
    easyExport(Cs_Table, "Cs", folder);
    easyExport(resis_Table, "resis", folder);
    easyExport(thermCon_Table, "thermCon", folder);

    easyExport(mass_frac_e_Table, "mass_frac_e", folder);
    easyExport(mass_frac_n_Table, "mass_frac_n", folder);
    easyExport(mass_frac_i_Table, "mass_frac_i", folder);

    cout << "EoS Constructed" << endl; 

    return 0; 
}
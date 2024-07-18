#include "mutation++.h"
#include "simulation.hpp"
#include "EoSUtils.hpp"
#include "parse.hpp"

using namespace Mutation;
using namespace Mutation::Thermodynamics;
using namespace Mutation::Transport; 
 
const int cellVarsNums = 12;

// SETTINGS ------

char delimiter = ' ';
double particleMass = protonMass;
string mixName = "HFusion";
string outputFolder = "../Calc/PIP/collRateStudy/";

// ROUNDING ------
double roundFactor = 1.0e8;

int nPoints = 5000;
int fixed_p = pAtmos;
int fixed_T = 8000;

int tempStart = 5000, tempEnd = 60000;
int pStart = 101325, pEnd = 10 * 101325;

string change_var = "p";

int main(void){

    // creating temp range ======

    vector<double> tempRange = linspace(tempStart, tempEnd, nPoints);
    vector<double> pRange = linspace(pStart, pEnd, nPoints);
    Mixture mix(mixName);

    if (change_var == "T"){
        outputFolder = outputFolder + "Tvar/";
    } else if (change_var == "p"){
        outputFolder = outputFolder + "pVar/";
    }

    // creating storage ======

    vector<double> mfn_array, mfi_array;
    vector<double> n_e_array, n_n_array, n_i_array;
    vector<double> coll_freq_in_array, coll_freq_en_array, coll_freq_ei_array, cLog_array;
    vector<double> alpha_in_array, alpha_en_array, alpha_ei_array;

    double temp, p;
    
    for (int i = 0; i < nPoints; i++){

        if (change_var == "T"){
            temp = tempRange[i];
            p = fixed_p;
        } else if (change_var == "p"){
            temp = fixed_T;
            p = pRange[i];
        }

        // EQUILIBRATE SYSTEM AND GET VALUES FROM Mutation++ ------
        mix.equilibrate(temp, p);

        double mass_frac_n = round(mix.Y()[1] * roundFactor) / roundFactor;
        double mass_frac_i = round(mix.Y()[2] * roundFactor) / roundFactor;

        double n_e = mix.numberDensity() * mix.X()[0];
        double n_n = mix.numberDensity() * mix.X()[1];
        double n_i = mix.numberDensity() * mix.X()[2];

        double cLog = getCoulombLog(n_e, temp);

        double coll_freq_in = get_coll_freq_in_noScale(n_n, particleMass, particleMass, temp);
        double coll_freq_en = get_coll_freq_en_noScale(n_n, particleMass, particleMass, temp);
        double coll_freq_ei = get_coll_freq_ei_noScale(n_i, cLog, temp);

        double alpha_in = get_coll_coeff_in_noScale(n_i, n_n, particleMass, particleMass, temp);
        double alpha_en = get_coll_coeff_en_noScale(n_i, n_n, particleMass, particleMass, temp);
        double alpha_ei = get_coll_coeff_ei_noScale(n_e, n_i, cLog, temp);

        mfi_array.push_back(mass_frac_i);
        mfn_array.push_back(mass_frac_n);
        n_e_array.push_back(n_e);
        n_n_array.push_back(n_n);
        n_i_array.push_back(n_i);
        cLog_array.push_back(cLog);
        coll_freq_ei_array.push_back(coll_freq_ei);
        coll_freq_en_array.push_back(coll_freq_en);
        coll_freq_in_array.push_back(coll_freq_in);
        alpha_in_array.push_back(alpha_in);
        alpha_en_array.push_back(alpha_en);
        alpha_ei_array.push_back(alpha_ei);

    }

    cout << "Exporting Results" << endl; 
    singleColExport(pRange, "pRange", outputFolder);
    singleColExport(tempRange, "tempRange", outputFolder);
    singleColExport(mfi_array, "mfi", outputFolder);
    singleColExport(mfn_array, "mfn", outputFolder);
    singleColExport(n_e_array, "n_e", outputFolder);
    singleColExport(n_i_array, "n_i", outputFolder);
    singleColExport(n_n_array, "n_n", outputFolder);
    singleColExport(cLog_array, "cLog", outputFolder);
    singleColExport(coll_freq_ei_array, "coll_freq_ei", outputFolder);
    singleColExport(coll_freq_en_array, "coll_freq_en", outputFolder);
    singleColExport(coll_freq_in_array, "coll_freq_in", outputFolder);
    singleColExport(alpha_in_array, "alpha_in", outputFolder);
    singleColExport(alpha_en_array, "alpha_en", outputFolder);
    singleColExport(alpha_ei_array, "alpha_ei", outputFolder);

    return 0; 
}
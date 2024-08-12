#include "mutation++.h"
#include "simulation.hpp"
#include "EoSUtils.hpp"
#include "parse.hpp"

using namespace Mutation;
using namespace Mutation::Thermodynamics;
using namespace Mutation::Transport; 
 
const int cellVarsNums = 12;
const int omp_threads = 4;

// SETTINGS ------

char delimiter = ' ';
double particleMass = protonMass;
string mixName = "HFusion";
string outputFolder = "../Calc/PIP/collRateStudy/p_n_T/";

// ROUNDING ------
double roundFactor = 1.0e8;

int nPoints_T = 1000;
int nPoints_P = 1000;

int tempStart = 5000, tempEnd = 20000;
int pStart = 1, pEnd = 0.1*101325;

string change_var = "p";

int main(void){

    // creating temp range ======

    vector<double> tempRange = linspace(tempStart, tempEnd, nPoints_T);
    vector<double> pRange = linspaceLog(pStart, pEnd, nPoints_P);
    Mixture mix(mixName);



    // creating storage ======

    Scalar2D mfi_table = makeScalar2D(nPoints_P, nPoints_T);
    Scalar2D mfn_table = makeScalar2D(nPoints_P, nPoints_T);
    Scalar2D n_e_table = makeScalar2D(nPoints_P, nPoints_T);
    Scalar2D n_n_table = makeScalar2D(nPoints_P, nPoints_T);
    Scalar2D n_i_table = makeScalar2D(nPoints_P, nPoints_T);
    Scalar2D coll_freq_in_table = makeScalar2D(nPoints_P, nPoints_T);
    Scalar2D coll_freq_en_table = makeScalar2D(nPoints_P, nPoints_T);
    Scalar2D coll_freq_ei_table = makeScalar2D(nPoints_P, nPoints_T);
    Scalar2D cLog_table = makeScalar2D(nPoints_P, nPoints_T);
    Scalar2D alpha_in_table = makeScalar2D(nPoints_P, nPoints_T);
    Scalar2D alpha_en_table = makeScalar2D(nPoints_P, nPoints_T);
    Scalar2D alpha_ei_table = makeScalar2D(nPoints_P, nPoints_T);

    Scalar2D gamma_table = makeScalar2D(nPoints_P, nPoints_T);
    Scalar2D resis_table = makeScalar2D(nPoints_P, nPoints_T);
    Scalar2D density_table = makeScalar2D(nPoints_P, nPoints_T);
    Scalar2D thermCon_ei_table = makeScalar2D(nPoints_P, nPoints_T);
    Scalar2D thermCon_in_table = makeScalar2D(nPoints_P, nPoints_T);



    double temp, p;
    
    for (int pIdx = 0; pIdx < nPoints_P; pIdx++){
        for (int tempIdx = 0; tempIdx < nPoints_T; tempIdx++){
            cout << "\rProgress: " << pIdx << "/" << nPoints_P << "     " <<std::flush;

            temp = tempRange[tempIdx];
            p = pRange[pIdx];

            // EQUILIBRATE SYSTEM AND GET VALUES FROM Mutation++ ------
            mix.equilibrate(temp, p);

            double mass_frac_n = round(mix.Y()[1] * roundFactor) / roundFactor;
            double mass_frac_i = round(mix.Y()[2] * roundFactor) / roundFactor;

            double density = mix.density();
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

            double gamma = mix.mixtureEquilibriumGamma();
            double resis = 1/mix.electricConductivity();

            double thermCon_ei = 3.9 * (n_i * kB * temp) / (coll_freq_ei * protonMass);
            double thermCon_in = 3.9 * (n_i * kB * temp) / (coll_freq_in * protonMass);


            mfi_table[pIdx][tempIdx] = mass_frac_i;
            mfn_table[pIdx][tempIdx] = mass_frac_n;
            n_e_table[pIdx][tempIdx] = n_e;
            n_n_table[pIdx][tempIdx] = n_n;
            n_i_table[pIdx][tempIdx] = n_i;
            cLog_table[pIdx][tempIdx] = cLog;
            coll_freq_ei_table[pIdx][tempIdx] = coll_freq_ei;
            coll_freq_en_table[pIdx][tempIdx] = coll_freq_en;
            coll_freq_in_table[pIdx][tempIdx] = coll_freq_in;
            alpha_in_table[pIdx][tempIdx] = alpha_in;
            alpha_en_table[pIdx][tempIdx] = alpha_en;
            alpha_ei_table[pIdx][tempIdx] = alpha_ei;
            gamma_table[pIdx][tempIdx] = gamma;
            resis_table[pIdx][tempIdx] = resis;
            density_table[pIdx][tempIdx] = density;
            thermCon_ei_table[pIdx][tempIdx] = thermCon_ei;
            thermCon_in_table[pIdx][tempIdx] = thermCon_in;

        }
    }

    cout << "Exporting Results" << endl; 
    singleColExport(pRange, "pRange", outputFolder);
    singleColExport(tempRange, "tempRange", outputFolder);
    easyExport(mfi_table, "mfi", outputFolder);
    easyExport(mfn_table, "mfn", outputFolder);
    easyExport(n_e_table, "n_e", outputFolder);
    easyExport(n_i_table, "n_i", outputFolder);
    easyExport(n_n_table, "n_n", outputFolder);
    easyExport(cLog_table, "cLog", outputFolder);
    easyExport(coll_freq_ei_table, "coll_freq_ei", outputFolder);
    easyExport(coll_freq_en_table, "coll_freq_en", outputFolder);
    easyExport(coll_freq_in_table, "coll_freq_in", outputFolder);
    easyExport(alpha_in_table, "alpha_in", outputFolder);
    easyExport(alpha_en_table, "alpha_en", outputFolder);
    easyExport(alpha_ei_table, "alpha_ei", outputFolder);
    easyExport(gamma_table, "gamma", outputFolder);
    easyExport(resis_table, "resis", outputFolder);
    easyExport(density_table, "density", outputFolder);
    easyExport(thermCon_ei_table, "thermCon_ei", outputFolder);
    easyExport(thermCon_in_table, "thermCon_in", outputFolder);
    

    return 0; 
}
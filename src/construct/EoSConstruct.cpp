#include "mutation++.h"
#include "simulation.hpp"
#include "EoSUtils.hpp"
#include "parse.hpp"

using namespace Mutation;
using namespace Mutation::Thermodynamics;
using namespace Mutation::Transport; 
 
const int cellVarsNums = 9;

// SETTINGS ------

char delimiter = ' ';
string folder = "./EoSData/HFusion4_1/";
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

    cout << "Constructing EoS " << endl;
    cout << "pLen: " << pLen << ", rhoLen: " << rhoLen << endl;

    // CREATING STORAGE FOR ALL VARIABLES ======
    Mixture mix("HFusion");

    Scalar2D pIndex = makeScalar2D(pLen, rhoLen);
    Scalar2D rhoIndex = makeScalar2D(pLen, rhoLen);

    Scalar2D T_Table = makeScalar2D(pLen, rhoLen);
    Scalar2D e_Table = makeScalar2D(pLen, rhoLen);
    Scalar2D Cs_Table = makeScalar2D(pLen, rhoLen);
    Scalar2D resis_Table = makeScalar2D(pLen, rhoLen);
    Scalar2D thermCon_Table = makeScalar2D(pLen, rhoLen);  

    Scalar2D n_e_Table = makeScalar2D(pLen, rhoLen);  
    Scalar2D n_n_Table = makeScalar2D(pLen, rhoLen);   
    Scalar2D n_i_Table = makeScalar2D(pLen, rhoLen);   
 
    // CONSTRUCTING ======
    for (int pIdx = 0; pIdx < pressures.size(); pIdx++){
        for (int rhoIdx = 0; rhoIdx < densities.size(); rhoIdx++){

            double pressure = pressures[pIdx];
            double rho = densities[rhoIdx];  
            double temp = BisectSolverEoS(rho, pressure, particleMass, mix, 1, 0.001, 200, 0.5);
            
            // cout << "temp: " << temp << endl;

            pIndex[pIdx][rhoIdx] = pressure * pScale;
            rhoIndex[pIdx][rhoIdx] = rho;
            T_Table[pIdx][rhoIdx] = temp;

            // EQUILIBRATE SYSTEM AND GET VALUES FROM Mutation++ ------
            mix.equilibrate(temp, pressure);

            /* Coefficients -> needs to be scaled */ 

            e_Table[pIdx][rhoIdx] = mix.mixtureEnergyMass() * eScale;
            Cs_Table[pIdx][rhoIdx] = mix.equilibriumSoundSpeed() * CsScale;
            resis_Table[pIdx][rhoIdx] = round(1/mix.electricConductivity() * resisScale * roundFactor)/roundFactor;
            thermCon_Table[pIdx][rhoIdx] = mix.equilibriumThermalConductivity() * thermConScale;

            /* Number densities */
            int numDensity = mix.numberDensity();
            n_e_Table[pIdx][rhoIdx] = round(numDensity * mix.X()[0] * roundFactor) / roundFactor;
            n_n_Table[pIdx][rhoIdx] = round(numDensity * mix.X()[1] * roundFactor) / roundFactor;
            n_i_Table[pIdx][rhoIdx] = round(numDensity * mix.X()[2] * roundFactor) / roundFactor;
        }
    }

    cout << "Exporting EoS" << endl; 

    easyExport(pIndex, "pressure", folder);
    easyExport(rhoIndex, "densities", folder);

    easyExport(T_Table, "T", folder);
    easyExport(e_Table, "e", folder);
    easyExport(Cs_Table, "Cs", folder);
    easyExport(resis_Table, "resis", folder);
    easyExport(thermCon_Table, "thermCon", folder);

    easyExport(n_e_Table, "n_e", folder);
    easyExport(n_n_Table, "n_n", folder);
    easyExport(n_i_Table, "n_i", folder);

    cout << "EoS Constructed" << endl; 

    return 0; 
}
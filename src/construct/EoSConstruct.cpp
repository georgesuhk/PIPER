#include "mutation++.h"
#include "simulation.hpp"
#include "parse.hpp"
#include "calcs.hpp"

using namespace Mutation;
using namespace Mutation::Thermodynamics;
using namespace Mutation::Transport; 
using namespace std;
 
const int cellVarsNums = 9;

// SETTINGS ------

char delimiter = ' ';
string folder = "../../EoSData/HFusion4_S1/";
double particleMass = protonMass;
string mixName = "HFusion";


int main(void){

    // LOADING IN P AND RHO RANGES ======

    vector<double> pressures = singleColParse(folder + "pRange.csv");
    vector<double> densities = singleColParse(folder + "rhoRange.csv");

    int pLen = pressures.size();
    int rhoLen = densities.size();

    cout << "Constructing EoS " << endl;
    cout << "pLen: " << pLen << ", rhoLen: " << rhoLen << endl;



    // CREATING STORAGE FOR ALL VARIABLES ======

    Mixture mix(mixName);

    Scalar2D pIndex = makeScalar2D(pLen, rhoLen);
    Scalar2D rhoIndex = makeScalar2D(pLen, rhoLen);

    Scalar2D T_Table = makeScalar2D(pLen, rhoLen);
    Scalar2D e_Table = makeScalar2D(pLen, rhoLen);
    Scalar2D Cs_Table = makeScalar2D(pLen, rhoLen);
    Scalar2D elecCon_Table = makeScalar2D(pLen, rhoLen);
    Scalar2D thermCon_Table = makeScalar2D(pLen, rhoLen);  
    Scalar2D thermCon_Table = makeScalar2D(pLen, rhoLen); 

    Scalar2D n_e_Table = makeScalar2D(pLen, rhoLen);  
    Scalar2D n_n_Table = makeScalar2D(pLen, rhoLen);   
    Scalar2D n_i_Table = makeScalar2D(pLen, rhoLen);   



    // CONSTRUCTING ======

    for (int pIdx = 0; pIdx < pressures.size(); pIdx++){
        for (int rhoIdx = 0; rhoIdx < densities.size(); rhoIdx ++){

        double pressure = pressures[pIdx];
        double density = densities[rhoIdx];  
        double temp = getIdealGasTemp(density, pressure, particleMass);

        pIndex[pIdx][rhoIdx] = pressure / pAtmos;
        T_Table[pIdx][rhoIdx] = temp;

        // EQUILIBRATE SYSTEM AND GET VALUES FROM Mutation++ ------

        mix.equilibrate(temp, pressure);

        /* Coefficients -> needs to be scaled */ 
        rhoIndex[pIdx][rhoIdx] = mix.density();
        e_Table[pIdx][rhoIdx] = mix.mixtureEnergyMass();
        Cs_Table[pIdx][rhoIdx] = mix.equilibriumSoundSpeed();
        elecCon_Table[pIdx][rhoIdx] = mix.electricConductivity();
        thermCon_Table[pIdx][rhoIdx] = mix.equilibriumThermalConductivity();

        /* Number densities */
        double numDensity = mix.numberDensity();
        n_e_Table[pIdx][rhoIdx] = numDensity * mix.X()[0];
        n_n_Table[pIdx][rhoIdx] = numDensity * mix.X()[1];
        n_i_Table[pIdx][rhoIdx] = numDensity * mix.X()[2];

        }
    }

    cout << "Exporting EoS" << endl; 

    easy2DExport(TVals, "temperature", outputFolder);
    easy2DExport(PVals, "pressure", outputFolder);

    easy2DExport(rhoMPP, "density", outputFolder);
    easy2DExport(eMPP, "e", outputFolder);
    easy2DExport(CsMPP, "Cs", outputFolder);
    easy2DExport(elecConMPP, "elecCon", outputFolder);
    easy2DExport(thermConMPP, "thermCon", outputFolder);
    easy2DExport(chargedMolFracMPP, "chargedMolFrac", outputFolder);

    return 0; 
}
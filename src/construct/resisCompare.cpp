#include "calcs.hpp"
#include "parse.hpp"
#include "mutation++.h"

using namespace Mutation;
using namespace Mutation::Thermodynamics;
using namespace Mutation::Transport; 
using namespace std;

int nPoints = 10000;
int p = 101325;
int tempStart = 5000, tempEnd = 60000;
// int tempStart = 30000, tempEnd = 300000;
const int cellVarsNums = 9;
Mixture mix("HFusion");

double m_e = electronMass;
double m_i = protonMass;
double m_n = protonMass;

string outputFolder = "../Calc/PIP/resisStudy/";
 
int main(void){

    int nSpecies = mix.nSpecies();

    vector<double> tempRange = linspace(tempStart, tempEnd, nPoints);
    vector<double> resis_MPP1_array, resis_MPP2_array, resis_MPP3_array, resis_brag_array;
    vector<double> n_e_array, n_i_array, n_n_array;
    vector<double> alpha_ei_array, alpha_en_array, alpha_in_array;
    vector<double> coulombLogArray;
    
    double T, numDensity;
    double n_e, n_i, n_n, alpha_ei, alpha_en, alpha_in, coulombLog;
    double resisMPP1, resisMPP2, resisMPP3, resisBrag;

    for (int i = 0; i < nPoints; i++){
        T = tempRange[i];

        mix.equilibrate(T, p);
        numDensity = mix.numberDensity();

        // resis from MPP ------
        resisMPP1 = 1 / mix.electricConductivity(1);
        resisMPP2 = 1 / mix.electricConductivity(2);
        resisMPP3 = 1 / mix.electricConductivity(3);


        // resis from Braginskii / general theory ------ 

        // species order e-, n, i
        n_e = numDensity * mix.X()[0];
        n_n = numDensity * mix.X()[1];
        n_i = numDensity * mix.X()[2];

        coulombLog = getCoulombLog(n_e, T);

        alpha_ei = get_coll_coeff_ei(n_e, n_i, coulombLog, T);
        alpha_en = get_coll_coeff_en(n_i, n_n, m_i, m_n, T);
        alpha_in = get_coll_coeff_in(n_i, n_n, m_i, m_n, T);

        
        resisBrag = (alpha_ei + alpha_en) / (pow(elementaryCharge, 2.0) * pow((n_e), 2.0));

        // STORING ======
        coulombLogArray.push_back(coulombLog);

        n_e_array.push_back(n_e);
        n_n_array.push_back(n_n);
        n_i_array.push_back(n_i);

        alpha_ei_array.push_back(alpha_ei);
        alpha_en_array.push_back(alpha_en);
        alpha_in_array.push_back(alpha_in);

        resis_MPP1_array.push_back(resisMPP1);
        resis_MPP2_array.push_back(resisMPP2);
        resis_MPP3_array.push_back(resisMPP3);

        resis_brag_array.push_back(resisBrag);
    }

    singleColExport(tempRange, "tempRange", outputFolder);

    singleColExport(coulombLogArray, "cLog", outputFolder);

    singleColExport(n_e_array, "n_e", outputFolder);
    singleColExport(n_n_array, "n_n", outputFolder);
    singleColExport(n_i_array, "n_i", outputFolder);

    singleColExport(alpha_ei_array, "alpha_ei", outputFolder);
    singleColExport(alpha_en_array, "alpha_en", outputFolder);
    singleColExport(alpha_in_array, "alpha_in", outputFolder);

    singleColExport(resis_MPP1_array, "resis_MPP1", outputFolder);
    singleColExport(resis_MPP2_array, "resis_MPP2", outputFolder);
    singleColExport(resis_MPP3_array, "resis_MPP3", outputFolder);
    singleColExport(resis_brag_array, "resis_Brag", outputFolder);






    return 0; 
}

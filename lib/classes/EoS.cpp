#include "EoS.hpp"

// BASE EOS CLASS ======

string EoS::getEoSType(){
    return EoSType;
}

double EoS::interp_Resis(double& rho, double& p){
    throw runtime_error("interp_Resis() should only be run if EoS is of type 'TabEoS'");
};




// IDEAL GAS EOS ======

double IdealEoS::get_gamma(){
    return gamma;
}

double IdealEoS::get_e(double& rho, double& p){
    return (p / (gamma - 1)) / rho;
};

double IdealEoS::get_T(double& rho, double& p){
    double n = rho / m;
    return p / (n * kB);
};

double IdealEoS::get_p(double& rho, double& e, bool interp){
    return (gamma - 1) * rho * e;
};

double IdealEoS::get_Cs(double& rho, double& p){
    return sqrt(gamma * p / rho);
};

double IdealEoS::get_Resis(int& i, int& j){
    return constResis;
};

double IdealEoS::interp_Resis(double& rho, double& p){
    return constResis;
}

void IdealEoS::set_gamma(double inputGamma){
    gamma = inputGamma;
};

void IdealEoS::set_constResis(double inputResis){
    constResis = inputResis;
};




// TAB EOS ======

// Internal solvers ------

int TabEoS::get_pLen(){
    return pressures.size();
}

int TabEoS::get_rhoLen(){
    return densities.size();
}

double TabEoS::get_p(double& rho, double& e, bool interp){
    int rho_lower_idx = getLowerBound(rho, activeRhoIndices, densities);
    return BisectSolver(rho, rho_lower_idx, e, densities, e_Table, pressures, 1e-5, 200);  
};

double TabEoS::get_e(double& rho, double& p){
    int p_lower_idx = getLowerBound(p, activePIndices, pressures);
    int rho_lower_idx = getLowerBound(rho, activeRhoIndices, densities);

    double p_lower = pressures[p_lower_idx];
    double p_higher = pressures[p_lower_idx + 1];

    double rho_lower = densities[rho_lower_idx];
    double rho_higher = densities[rho_lower_idx + 1];

    return bilinearInterp(e_Table, rho, p, rho_lower, rho_higher, rho_lower_idx, p_lower, p_higher, p_lower_idx);
}

double TabEoS::get_T(double& rho, double& p){
    int p_lower_idx = getLowerBound(p, activePIndices, pressures);
    int rho_lower_idx = getLowerBound(rho, activeRhoIndices, densities);

    double p_lower = pressures[p_lower_idx];
    double p_higher = pressures[p_lower_idx + 1];

    double rho_lower = densities[rho_lower_idx];
    double rho_higher = densities[rho_lower_idx + 1];

    return bilinearInterp(T_Table, rho, p, rho_lower, rho_higher, rho_lower_idx, p_lower, p_higher, p_lower_idx);
}

double TabEoS::get_Cs(double& rho, double& p){
    int p_lower_idx = getLowerBound(p, activePIndices, pressures);
    int rho_lower_idx = getLowerBound(rho, activeRhoIndices, densities);

    double p_lower = pressures[p_lower_idx];
    double p_higher = pressures[p_lower_idx + 1];

    double rho_lower = densities[rho_lower_idx];
    double rho_higher = densities[rho_lower_idx + 1];

    return bilinearInterp(Cs_Table, rho, p, rho_lower, rho_higher, rho_lower_idx, p_lower, p_higher, p_lower_idx);
}

double TabEoS::get_Resis(int& i, int& j){

    // will also need to apply scaling here
    return 0 * resisScaling;
}

double TabEoS::interp_Resis(double& rho, double& p){
    int p_lower_idx = getLowerBound(p, activePIndices, pressures);
    int rho_lower_idx = getLowerBound(rho, activeRhoIndices, densities);

    double p_lower = pressures[p_lower_idx];
    double p_higher = pressures[p_lower_idx + 1];

    double rho_lower = densities[rho_lower_idx];
    double rho_higher = densities[rho_lower_idx + 1];

    double elecConduct = bilinearInterp(elecConduct_Table, rho, p, rho_lower, rho_higher, rho_lower_idx, p_lower, p_higher, p_lower_idx);
    return 1/elecConduct * resisScaling;

}

void TabEoS::genFromData(Mesh2D mesh, vector<string> varList, string dataFolder, char delimiter){
    cout << "Loading in TabEoS data from: " << dataFolder << endl;

    Scalar2D tabularData;
    resis_Cache = makeScalar2D(mesh.nCellsX+2, mesh.nCellsY+2);

    for (string& var : varList){
        string filename = dataFolder + var + ".csv";
        tabularData = table_from_csv(filename, delimiter);

        if (var == "pressure"){
            pressures = getColumnEoS(tabularData, 1);
        }
        else if (var == "density"){
            densities = tabularData[0];
        }
        else if (var == "temperature"){
            T_Table = tabularData;
        }
        else if (var == "Cs"){
            Cs_Table = tabularData;
        }
        else if (var == "e"){
            e_Table = tabularData;
        }
        else if (var == "elecCon"){
            elecConduct_Table = tabularData;
        }
        else if (var == "thermCon"){
            thermConduct_Table = tabularData;
        }
        else {
            cout << "Warning: Variable " << var << " is invalid.";
        }
    }

    activeRhoIndices = {0, int(densities.size())-1};
    activePIndices = {0, int(pressures.size())-1};

    cout << "TabEoS Loaded. \n" << endl;
}
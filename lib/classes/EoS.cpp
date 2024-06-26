#include "EoS.hpp"

// BASE EOS CLASS ======

string EoS::getEoSType(){
    return EoSType;
}

double EoS::interp_Resis(double& rho, double& p){
    throw runtime_error("interp_Resis() should only be run if EoS is of type 'TabEoS'");
};

void EoS::cache_p(Vec2D& u, Mesh2D& mesh){
    if (p_cache.empty()){
        p_cache = makeScalar2D(mesh.nCellsX+2, mesh.nCellsY+2);
    }

    double rho, KE, MagE, e;
    for (int i = 1; i < mesh.nCellsX+1; i++){
        for (int j = 1; j < mesh.nCellsY+1; j++){
            rho = u[i][j][0];
            KE = 0.5 * getVectorMagnitude({u[i][j][1]/rho, u[i][j][2]/rho, u[i][j][3]/rho});
            MagE = 0.5 * getVectorMagnitude({u[i][j][5], u[i][j][6], u[i][j][7]});
            
            e = (u[i][j][4] - KE - MagE)/rho;
            p_cache[i][j] = interp_p(rho, e);
        }
    }
}

void EoS::cache_resis(Vec2D& u, Mesh2D& mesh){
    if (resis_cache.empty()){
        cout << "empty" << endl;
        resis_cache = makeScalar2D(mesh.nCellsX+2, mesh.nCellsY+2);
    }

    double rho, p;
    for (int i = 1; i < mesh.nCellsX+1; i++){
        for (int j = 1; j < mesh.nCellsY+1; j++){
            rho = u[i][j][0];
            p = get_p(i,j);
            
            resis_cache[i][j] = interp_Resis(rho, p);
        }
    }
}

/* current caches p and resis */
void EoS::cacheAll(Vec2D& u, Mesh2D& mesh){
    if (p_cache.empty()){
        // contains ghost cells as usual
        p_cache = makeScalar2D(mesh.nCellsX+2, mesh.nCellsY+2);
    }

    if (resis_cache.empty()){
        // contains ghost cells as usual
        resis_cache = makeScalar2D(mesh.nCellsX+2, mesh.nCellsY+2);
    }

    double rho, vx, vy, vz, e, KE, BMag, MagE, p;
    for (int i = 1; i < mesh.nCellsX+1; i++){
        for (int j = 1; j < mesh.nCellsY+1; j++){
            rho = u[i][j][0];
            vx = u[i][j][1] / rho;
            vy = u[i][j][2] / rho;
            vz = u[i][j][3] / rho;

            KE = 0.5 * rho * ( pow(vx, 2.0) + pow(vy, 2.0) + pow(vz, 2.0) );
            BMag = getVectorMagnitude({u[i][j][5], u[i][j][6], u[i][j][7]});
            MagE = 0.5 * pow(BMag, 2.0);
            e = (u[i][j][4] - KE - MagE)/rho;
            p = interp_p(rho, e);
            
            resis_cache[i][j] = interp_Resis(rho, p);
            p_cache[i][j] = p;
        }
    }
}

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

double IdealEoS::get_p(int& i, int& j){
    return p_cache[i][j];
};

double IdealEoS::interp_p(double& rho, double& e){
    return (gamma - 1) * rho * e;
}

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

Scalar1D& TabEoS::get_pressures(){
    return pressures;
}

Scalar1D& TabEoS::get_densities(){
    return densities;
}

Scalar2D& TabEoS::get_T_Table(){
    return T_Table;
}

int TabEoS::get_pLen(){
    return pressures.size();
}

int TabEoS::get_rhoLen(){
    return densities.size();
}

double TabEoS::get_p(int& i, int& j){
    return p_cache[i][j];
}

double TabEoS::interp_p(double& rho, double& e){
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
    return resis_cache[i][j];
}

double TabEoS::interp_Resis(double& rho, double& p){
    int p_lower_idx = getLowerBound(p, activePIndices, pressures);
    int rho_lower_idx = getLowerBound(rho, activeRhoIndices, densities);

    double p_lower = pressures[p_lower_idx];
    double p_higher = pressures[p_lower_idx + 1];

    double rho_lower = densities[rho_lower_idx];
    double rho_higher = densities[rho_lower_idx + 1];

    double resis = bilinearInterp(resis_Table, rho, p, rho_lower, rho_higher, rho_lower_idx, p_lower, p_higher, p_lower_idx);
    return resis;

}

void TabEoS::genFromData(Mesh2D mesh, vector<string> varList, string dataFolder, char delimiter){
    cout << "Loading in TabEoS data from: " << dataFolder << endl;

    Scalar2D tabularData;

    for (string& var : varList){
        string filename = dataFolder + var + ".csv";
        tabularData = easyParseTable(filename, delimiter);

        if (var == "pressure"){
            pressures = getColumnEoS(tabularData, 1);
        }
        else if (var == "densities"){
            densities = tabularData[800];
        }
        else if (var == "T"){
            T_Table = tabularData;
        }
        else if (var == "Cs"){
            Cs_Table = tabularData;
        }
        else if (var == "e"){
            e_Table = tabularData;
        }
        else if (var == "resis"){
            resis_Table = tabularData;
        }
        else if (var == "thermCon"){
            thermConduct_Table = tabularData;
        }
        else if (var == "n_e"){
            n_e_Table = tabularData;
        }
        else if (var == "n_n"){
            n_n_Table = tabularData;
        }
        else if (var == "n_i"){
            n_i_Table = tabularData;
        }
        else {
            cout << "Warning: Variable " << var << " is invalid." << endl;
        }
    }

    activeRhoIndices = {0, int(densities.size())-1};
    activePIndices = {0, int(pressures.size())-1};

    cout << "rho last: " << densities.back() << endl;
    cout << "TabEoS Loaded. \n" << endl;
}
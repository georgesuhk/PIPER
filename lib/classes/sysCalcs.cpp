#include "sysCalcs.hpp"

// SYSCALCS (Base) FUNCTIONS ======

// operated on CellVec in primitive form ------
double SysCalcs::get_e_Prim(CellVec& uPrim){
    double rho = uPrim[0];
    double p = uPrim[4];
    return EoSPtr->get_e(rho, p);
};

double SysCalcs::get_T_Prim(CellVec& uPrim){

    return EoSPtr->get_T(uPrim[0], uPrim[4]);
}

double SysCalcs::get_KE_Prim(CellVec& uPrim){

    return 0.5 * uPrim[0] * ( pow(uPrim[1], 2.0) + pow(uPrim[2], 2.0) + pow(uPrim[3], 2.0) );
};

double SysCalcs::get_E_Prim(CellVec& uPrim){

    double rho = uPrim[0];
    double p = uPrim[4];
    double e = EoSPtr->get_e(rho, p);
    double KE = get_KE_Prim(uPrim);
    return rho * e + KE;  
}

double SysCalcs::get_U(CellVec& uPrim){
    return get_E_Prim(uPrim) + get_MagE(uPrim);
};

double SysCalcs::interp_mass_frac_n_Prim(CellVec& uPrim){
    double rho = uPrim[0];
    double p = uPrim[4];
    return EoSPtr->interp_mass_frac_n(rho, p);
}

double SysCalcs::interp_mass_frac_i_Prim(CellVec& uPrim){
    double rho = uPrim[0];
    double p = uPrim[4];
    return EoSPtr->interp_mass_frac_i(rho, p);
}


// operated on CellVec in conservative form ------

double SysCalcs::get_T(CellVec& u){
    double p = get_p(u);
    return EoSPtr->get_T(u[0], p);
}

double SysCalcs::get_KE(CellVec& u){

    double rho = u[0];
    double vx = u[1] / rho;
    double vy = u[2] / rho;
    double vz = u[3] / rho;

    return 0.5 * rho * ( pow(vx, 2.0) + pow(vy, 2.0) + pow(vz, 2.0) );
};

double SysCalcs::get_KE(double& rho, double& vx, double& vy, double& vz){
    return 0.5 * rho * ( pow(vx, 2.0) + pow(vy, 2.0) + pow(vz, 2.0) );   
};

double SysCalcs::get_E(CellVec& u){
    return u[4] - get_MagE(u);  
}


double SysCalcs::get_e(CellVec& u){
    return (get_E(u) - get_KE(u)) / u[0];
};

double SysCalcs::interp_gamma(CellVec& u){
    double rho = u[0];
    double p = get_p(u);
    return EoSPtr->get_gamma(rho, p);
};

double SysCalcs::get_p(CellVec& u){
    double e = get_e(u);
    return EoSPtr->interp_p(u[0], e);
};

double SysCalcs::get_p(int i, int j){
    return EoSPtr->get_p(i, j);
}

double SysCalcs::get_m_i(){
    return m_i;
}

double SysCalcs::get_m_n(){
    return m_n;
}

double SysCalcs::interp_mass_frac_e(CellVec& u){
    double rho = u[0];
    double p = get_p(u);
    return EoSPtr->interp_mass_frac_e(rho, p);
}

double SysCalcs::interp_mass_frac_n(CellVec& u){
    double rho = u[0];
    double p = get_p(u);
    double mass_frac_n = EoSPtr->interp_mass_frac_n(rho, p);
    return mass_frac_n;
}

double SysCalcs::interp_mass_frac_i(CellVec& u){
    double rho = u[0];
    double p = get_p(u);
    double mass_frac_i = EoSPtr->interp_mass_frac_i(rho, p);
    if (mass_frac_i == 0){
        mass_frac_i += 1e-6;
    }
    return mass_frac_i;
}

double SysCalcs::get_n_n(double rho, double mass_frac_n){
    return mass_frac_n * rho / m_n;
}

double SysCalcs::get_n_i(double rho, double mass_frac_i){
    return mass_frac_i * rho / m_i;
}

double SysCalcs::get_Resis(CellVec& u, int i, int j, bool interp){
    double resis;
    double rho = u[0];
    double p = get_p(u);

    if (interp){
        // re-interprelating the resistivity for current cell
        resis = EoSPtr->interp_Resis(rho, p);
        
    } else {
        // accessing cached variable instead
        resis = EoSPtr->get_Resis(i, j); 
    }

    return resis;
}


// operated on both prim and conserv ------

double SysCalcs::get_MagE(CellVec& u){

    double BMag = getVectorMagnitude({u[5], u[6], u[7]});
    return 0.5 * pow(BMag, 2.0);
};



// Wave speed estimations ======
double SysCalcs::get_Cs(CellVec& u){
    double p = get_p(u);
    return EoSPtr->get_Cs(u[0], p);
};

double SysCalcs::get_Ca(CellVec& u, char axis){

    double Bk;
    switch (axis){
        case 'x':
            Bk = u[5];
            break;
        case 'y':
            Bk = u[6];
            break;
        case 'z':
            Bk = u[7];
            break;

    }

    return fabs(Bk)/sqrt(u[0]);    
};

double SysCalcs::get_Cf(CellVec& u, char axis){

    double Bk;
    switch (axis){
        case 'x':
            Bk = u[5];
            break;
        case 'y':
            Bk = u[6];
            break;
        case 'z':
            Bk = u[7];
            break;
        default:
            throw invalid_argument("Invalid axis encountered in getCf()");
    }

    double Cs = get_Cs(u);
    double Ca = get_Ca(u, axis);

    double innerSqrt;
    innerSqrt = sqrt( pow(pow(Cs, 2.0) + pow(Ca, 2.0), 2.0) - 4 * (pow(Cs, 2.0) * pow(Bk, 2.0)) / u[0] ); 

    return sqrt( 0.5 * ( pow(Cs, 2.0) + pow(Ca, 2.0) + innerSqrt ));  
};

string SysCalcs::getSysName(){
    return sysName;
}

shared_ptr<EoS> SysCalcs::getEoSPtr(){
    return EoSPtr;
}




// FULLY IONISED PLASMA SYSTEM CALCS ======

CellVec FIPCalcs::primToConserv(CellVec& uPrim){
    CellVec uConserv(cellVecLen, 0);

    uConserv[0] = uPrim[0];
    uConserv[1] = uPrim[0] * uPrim[1];
    uConserv[2] = uPrim[0] * uPrim[2];
    uConserv[3] = uPrim[0] * uPrim[3];
    uConserv[4] = get_U(uPrim);
    uConserv[5] = uPrim[5];
    uConserv[6] = uPrim[6];
    uConserv[7] = uPrim[7];
    uConserv[8] = uPrim[8];

    return uConserv;
};

CellVec FIPCalcs::conservToPrim(CellVec& u){
    CellVec uPrim(cellVecLen, 0);

    uPrim[0] = u[0];
    uPrim[1] = u[1] / u[0];
    uPrim[2] = u[2] / u[0];
    uPrim[3] = u[3] / u[0];
    uPrim[4] = get_p(u);
    uPrim[5] = u[5];
    uPrim[6] = u[6];
    uPrim[7] = u[7];
    uPrim[8] = u[8];

    return uPrim;
};

CellVec FIPCalcs::f(CellVec& u){
    CellVec f_out(cellVecLen, 0);

    double rho = u[0];
    double vx = u[1] / rho;
    double vy = u[2] / rho;
    double vz = u[3] / rho;
    double U = u[4];
    double Bx = u[5];
    double By = u[6];
    double Bz = u[7];
    double psi = u[8];
    
    double p = get_p(u);
    double MagE = get_MagE(u);
    double vDotB = dotProduct({vx, vy, vz}, {Bx, By, Bz});

    f_out[0] = rho * vx;
    f_out[1] = rho * pow(vx, 2.0) + p + MagE - pow(Bx, 2.0);
    f_out[2] = rho * vx * vy - Bx * By;
    f_out[3] = rho * vx * vz - Bx * Bz;
    f_out[4] = (U + p + MagE) * vx - vDotB * Bx;
    f_out[5] = psi;
    f_out[6] = By * vx - Bx * vy;
    f_out[7] = Bz * vx - Bx * vz;
    f_out[8] = 0;

    return f_out;
};

CellVec FIPCalcs::g(CellVec& u){
    CellVec g_out(cellVecLen, 0);

    double rho = u[0];
    double vx = u[1] / rho;
    double vy = u[2] / rho;
    double vz = u[3] / rho;
    double U = u[4];
    double Bx = u[5];
    double By = u[6];
    double Bz = u[7];
    double psi = u[8];

    double p = get_p(u);
    double MagE = get_MagE(u);
    double vDotB = dotProduct({vx, vy, vz}, {Bx, By, Bz});

    g_out[0] = rho * vy;
    g_out[1] = rho * vy * vx - By * Bx;
    g_out[2] = rho * pow(vy, 2.0) + p + MagE - pow(By, 2.0);
    g_out[3] = rho * vy * vz - By * Bz;
    g_out[4] = (U + p + MagE) * vy - vDotB * By;
    g_out[5] = Bx * vy - By * vx;
    g_out[6] = psi;
    g_out[7] = Bz * vy - By * vz;
    g_out[8] = 0;

    return g_out;
};

double FIPCalcs::getFastestWaveSpeed(CellVec& u){

    double vx = u[1] / u[0];
    double vy = u[2] / u[0];
    double vz = u[3] / u[0];

    double a_x = fabs(vx) + get_Cf(u, 'x');
    double a_y = fabs(vy) + get_Cf(u, 'y');
    double a_z = fabs(vz) + get_Cf(u, 'z');

    double a_array[] = {a_x, a_y, a_z}; 
    double a = *max_element(begin(a_array), end(a_array));

    return a;
};



// PIP0 CALCS ======


CellVec PIP0_Calcs::primToConserv(CellVec& uPrim){
    CellVec uConserv(cellVecLen, 0);

    uConserv[0] = uPrim[0];
    uConserv[1] = uPrim[0] * uPrim[1];
    uConserv[2] = uPrim[0] * uPrim[2];
    uConserv[3] = uPrim[0] * uPrim[3];
    uConserv[4] = get_U(uPrim);
    uConserv[5] = uPrim[5];
    uConserv[6] = uPrim[6];
    uConserv[7] = uPrim[7];
    uConserv[8] = uPrim[8];
    // w values
    uConserv[9] = uPrim[9]; // wx
    uConserv[10] = uPrim[10]; // wy
    uConserv[11] = uPrim[11]; // wz


    return uConserv;
};

double PIP0_Calcs::get_total_neutral_E(CellVec& u){
    double rho = u[0];
    double wx = u[9];
    double wy = u[10];
    double wz = u[11];
    double mass_frac_i = interp_mass_frac_i(u);
    double mass_frac_n = interp_mass_frac_n(u);
    double n_n = get_n_n(rho, mass_frac_n);
    double gamma = interp_gamma(u);
    double rho_n = rho * mass_frac_n;
    
    // getting v_n (neutral speed)
    double vn_x = u[1]/rho + mass_frac_i * wx;
    double vn_y = u[2]/rho + mass_frac_i * wy;
    double vn_z = u[3]/rho + mass_frac_i * wz;

    double KE_n = get_KE(rho_n, vn_x, vn_y, vn_z);
    double p_n = n_n * kBScaled * get_T(u);
    double e_n_plus_p_n = (gamma)/(gamma-1) * p_n;

    return KE_n + e_n_plus_p_n;
}

CellVec PIP0_Calcs::conservToPrim(CellVec& u){
    CellVec uPrim(cellVecLen, 0);

    uPrim[0] = u[0];
    uPrim[1] = u[1] / u[0];
    uPrim[2] = u[2] / u[0];
    uPrim[3] = u[3] / u[0];
    uPrim[4] = get_p(u);
    uPrim[5] = u[5];
    uPrim[6] = u[6];
    uPrim[7] = u[7];
    uPrim[8] = u[8];

    // w values
    uPrim[9] = u[9]; // wx
    uPrim[10] = u[10]; // wy
    uPrim[11] = u[11]; // wz

    return uPrim;
};

CellVec PIP0_Calcs::f(CellVec& u){
    CellVec f_out(cellVecLen, 0);

    double rho = u[0];
    double vx = u[1] / rho;
    double vy = u[2] / rho;
    double vz = u[3] / rho;
    double U = u[4];
    double Bx = u[5];
    double By = u[6];
    double Bz = u[7];
    double psi = u[8];
    double wx = u[9];
    double wy = u[10];
    double wz = u[11];
    
    double p = get_p(u);
    double E = get_E(u);
    double total_neutral_E = get_total_neutral_E(u);
    double MagE = get_MagE(u);
    double vDotB = dotProduct({vx, vy, vz}, {Bx, By, Bz});

    // PIP modifications
    double mass_frac_n = interp_mass_frac_n(u);
    double mass_frac_i = interp_mass_frac_i(u);

    f_out[0] = rho * vx;
    f_out[1] = rho * pow(vx, 2.0) - rho * mass_frac_i * mass_frac_n * pow(wx, 2.0) + p + MagE - pow(Bx, 2.0);
    f_out[2] = rho * vx * vy - rho * mass_frac_i * mass_frac_n * wx * wy - Bx * By;
    f_out[3] = rho * vx * vz - rho * mass_frac_i * mass_frac_n * wx * wz - Bx * Bz;
    f_out[4] = (U + p + MagE) * vx + (mass_frac_n * (E + p) - total_neutral_E) * wx - vDotB * Bx; // need to fix energy conservation
    f_out[5] = psi;
    f_out[6] = By * vx - Bx * vy;
    f_out[7] = Bz * vx - Bx * vz;
    f_out[8] = 0;
    // w flux
    // f_out[9] = 0;
    // f_out[10] = 0;
    // f_out[11] = 0;
    f_out[9] = vx * wx + mass_frac_n * pow(wx, 2.0) - (1 / (mass_frac_i * rho)) * (MagE - pow(Bx, 2.0));
    f_out[10] = vx * wy + mass_frac_n * wx * wy + (1 / (mass_frac_i * rho)) * (Bx * By);
    f_out[11] = vx * wz + mass_frac_n * wx * wz + (1 / (mass_frac_i * rho)) * (Bx * Bz);

    return f_out;
};

CellVec PIP0_Calcs::g(CellVec& u){
    CellVec g_out(cellVecLen, 0);

    double rho = u[0];
    double vx = u[1] / rho;
    double vy = u[2] / rho;
    double vz = u[3] / rho;
    double U = u[4];
    double Bx = u[5];
    double By = u[6];
    double Bz = u[7];
    double psi = u[8];
    double wx = u[9];
    double wy = u[10];
    double wz = u[11];

    double p = get_p(u);
    double E = get_E(u);
    double total_neutral_E = get_total_neutral_E(u);
    double MagE = get_MagE(u);
    double vDotB = dotProduct({vx, vy, vz}, {Bx, By, Bz});

    // PIP modifications
    double mass_frac_n = interp_mass_frac_n(u);
    double mass_frac_i = interp_mass_frac_i(u);

    g_out[0] = rho * vy;
    g_out[1] = rho * vy * vx - By * Bx - rho * mass_frac_i * mass_frac_n * wy * wx;
    g_out[2] = rho * pow(vy, 2.0) - rho * mass_frac_i * mass_frac_n * pow(wy, 2.0) + p + MagE - pow(By, 2.0);
    g_out[3] = rho * vy * vz - By * Bz - rho * mass_frac_i * mass_frac_n * wy * wz;
    g_out[4] = (U + p + MagE) * vy + (mass_frac_n * (E + p) - total_neutral_E) * wy - vDotB * By;
    g_out[5] = Bx * vy - By * vx;
    g_out[6] = psi;
    g_out[7] = Bz * vy - By * vz;
    g_out[8] = 0;
    // w flux
    // g_out[9] = 0;
    // g_out[10] = 0;
    // g_out[11] = 0;
    g_out[9] = vy * wx + mass_frac_n * wy * wx + (1 / (mass_frac_i * rho)) * (By * Bx);
    g_out[10] = vy * wy + mass_frac_n * pow(wy, 2.0) - (1 / (mass_frac_i * rho)) * (MagE - pow(By, 2.0));
    g_out[11] = vy * wz + mass_frac_n * wy * wz + (1 / (mass_frac_i * rho)) * (By * Bz);

    return g_out;
};

// fastest wave speed taking into account of w
double PIP0_Calcs::getFastestWaveSpeed(CellVec& u){

    double vx = u[1] / u[0];
    double vy = u[2] / u[0];
    double vz = u[3] / u[0];
    double wx = u[9];
    double wy = u[10];
    double wz = u[11];

    double a_x = fabs(vx) + get_Cf(u, 'x');
    double a_y = fabs(vy) + get_Cf(u, 'y');
    double a_z = fabs(vz) + get_Cf(u, 'z');
    double aw_x = fabs(wx) + fabs(a_x);
    double aw_y = fabs(wy) + fabs(a_y);
    double aw_z = fabs(wz) + fabs(a_z);

    double a_array[] = {aw_x, aw_y, aw_z}; 
    double a = *max_element(begin(a_array), end(a_array));

    return a;
};






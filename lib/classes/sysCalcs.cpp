#include "sysCalcs.hpp"

// SYSCALCS (Base) FUNCTIONS ======

// operated on CellVec in primitive form ------
double SysCalcs::get_e_Prim(CellVec& uPrim){
    double rho = uPrim[0];
    double p = uPrim[4];
    return EoSPtr->get_e(rho, p);
};

double SysCalcs::get_T_Prim(CellVec& uPrim){

    return EoSPtr->interp_T(uPrim[0], uPrim[4]);
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

double SysCalcs::get_T(int i, int j){
    return EoSPtr->get_T(i, j);
}

double SysCalcs::interp_T(CellVec& u){
    double rho = u[0];
    double p = interp_p(u);
    return EoSPtr->interp_T(rho, p);
}

double SysCalcs::interp_T(double& rho, double& p){
    return EoSPtr->interp_T(rho, p);
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
    double p = interp_p(u);
    return EoSPtr->interp_gamma(rho, p);
};

double SysCalcs::interp_p(CellVec& u){
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
    double p = interp_p(u);
    return EoSPtr->interp_mass_frac_e(rho, p);
}

double SysCalcs::interp_mass_frac_n(CellVec& u){
    double rho = u[0];
    double p = interp_p(u);
    double mass_frac_n = EoSPtr->interp_mass_frac_n(rho, p);
    if (mass_frac_n == 0){
        mass_frac_n += 1e-200;
    }
    return mass_frac_n;
}

double SysCalcs::get_mass_frac_i(int i, int j){
    return EoSPtr->get_mass_frac_i(i, j);
}

double SysCalcs::interp_mass_frac_i(CellVec& u){
    double rho = u[0];
    double p = interp_p(u);
    double mass_frac_i = EoSPtr->interp_mass_frac_i(rho, p);
    if (mass_frac_i == 0){
        mass_frac_i += 1e-200;
    }
    return mass_frac_i;
}

double SysCalcs::interp_mass_frac_i(double& rho, double& p){
    return EoSPtr->interp_mass_frac_i(rho, p);
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
    double p = interp_p(u);

    if (interp){
        // re-interprelating the resistivity for current cell
        resis = EoSPtr->interp_Resis(rho, p);
        
    } else {
        // accessing cached variable instead
        resis = EoSPtr->get_Resis(i, j); 
    }

    return resis;
}

double SysCalcs::interp_Resis(double& rho, double& p){
    return EoSPtr->interp_Resis(rho, p);
}

// operated on both prim and conserv ------

double SysCalcs::get_MagE(CellVec& u){

    double BMag = getVectorMagnitude({u[5], u[6], u[7]});
    return 0.5 * pow(BMag, 2.0);
};



// Wave speed estimations ======
double SysCalcs::get_Cs(CellVec& u){
    double p = interp_p(u);
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
    uPrim[4] = interp_p(u);
    uPrim[5] = u[5];
    uPrim[6] = u[6];
    uPrim[7] = u[7];
    uPrim[8] = u[8];

    return uPrim;
};

CellVec FIPCalcs::f(CellVec& u, int i, int j, bool interp){
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
    
    double p = interp_p(u);
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

CellVec FIPCalcs::g(CellVec& u, int i, int j, bool interp){
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

    double p = interp_p(u);
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
    double p_n = n_n * kBScaled * interp_T(u);
    double e_n_plus_p_n = (gamma)/(gamma-1) * p_n;

    return KE_n + e_n_plus_p_n;
}

CellVec PIP0_Calcs::conservToPrim(CellVec& u){
    CellVec uPrim(cellVecLen, 0);

    uPrim[0] = u[0];
    uPrim[1] = u[1] / u[0];
    uPrim[2] = u[2] / u[0];
    uPrim[3] = u[3] / u[0];
    uPrim[4] = interp_p(u);
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

CellVec PIP0_Calcs::f(CellVec& u, int i, int j, bool interp){
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
    
    double p = interp_p(u);
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
    f_out[9] = 0;
    f_out[10] = 0;
    f_out[11] = 0;


    // w evo with mag
    // f_out[9] = - (1 / (mass_frac_i * rho)) * (MagE - pow(Bx, 2.0));
    // f_out[10] = (1 / (mass_frac_i * rho)) * (Bx * By);
    // f_out[11] = (1 / (mass_frac_i * rho)) * (Bx * Bz);

    //w evo with no mag

    return f_out;
};

CellVec PIP0_Calcs::g(CellVec& u, int i, int j, bool interp){
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

    double p = interp_p(u);
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
    g_out[9] = 0;
    g_out[10] = 0;
    g_out[11] = 0;


    // g_out[9] = (1 / (mass_frac_i * rho)) * (By * Bx);
    // g_out[10] = - (1 / (mass_frac_i * rho)) * (MagE - pow(By, 2.0));
    // g_out[11] = (1 / (mass_frac_i * rho)) * (By * Bz);

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



void PIP0_Calcs::set_w_no_inert(Vec2D& u, Mesh2D& mesh){
    /* constant values */
    double two_dx = 2 * mesh.dx;
    double two_dy = 2 * mesh.dy;

    /* state vars */
    double Bx, By, Bz, mass_frac_i, mass_frac_n, n_i, n_n, m_i, m_n, T;

    /* state vars used in derivatives */
    double mfi_x_plus, mfi_x_minus, mfi_y_plus, mfi_y_minus;
    double T_x_plus, T_x_minus, T_y_plus, T_y_minus;

    /* derivatives */
    double dBx_dx, dBx_dy, dBy_dx, dBy_dy, dBz_dx, dBz_dy;
    double d_pn_dx, d_pn_dy, d_pei_dx, d_pei_dy;

    /* compound terms */
    double mag_prefac, mag_term_x, mag_term_y, mag_term_z;
    double J_term_prefac, J_x, J_y, J_z;
    double G_x, G_y, G_z;

    /* collision freqs */
    double alpha_n, alpha_en, alpha_in;

    for (int i = 1; i < mesh.nCellsX+1; i++){
        for (int j = 1; j < mesh.nCellsY+1; j++){

            /* setting state vars */
            Bx = u[i][j][5];
            By = u[i][j][6];
            Bz = u[i][j][7];

            mass_frac_i = interp_mass_frac_i(u[i][j]);
            mass_frac_n = 1 - mass_frac_i;
            n_i = get_n_i(u[i][j][0], mass_frac_i);
            n_n = get_n_n(u[i][j][0], mass_frac_n);
            m_n = get_m_n();
            m_i = get_m_i();
            T = get_T(i, j);

            /* calculating collision integrals */
            alpha_in = get_coll_coeff_in(n_i, n_n, m_i, m_n, T);
            alpha_en = get_coll_coeff_en(n_i, n_n, m_i, m_n, T);
            alpha_n = alpha_in + alpha_en;

            /* calculating derivatives */
            T_x_plus = get_T(i+1, j);
            T_x_minus = get_T(i-1, j);
            T_y_plus = get_T(i, j+1);
            T_y_minus = get_T(i, j-1);

            mass_frac_i = interp_mass_frac_i(u[i][j]);
            mass_frac_n = 1 - mass_frac_i;
            mfi_x_plus = interp_mass_frac_i(u[i+1][j]);
            mfi_x_minus = interp_mass_frac_i(u[i-1][j]);
            mfi_y_plus = interp_mass_frac_i(u[i][j+1]);
            mfi_y_minus = interp_mass_frac_i(u[i][j-1]);

            /* partial pressure derivatives */
            d_pn_dx = kBScaled * ( get_n_n(u[i+1][j][0], 1-mfi_x_plus) * T_x_plus - get_n_n(u[i-1][j][0], 1-mfi_x_minus) * T_x_minus) / two_dx;
            d_pn_dy = kBScaled * ( get_n_n(u[i][j+1][0], 1-mfi_y_plus) * T_y_plus - get_n_n(u[i][j-1][0], 1-mfi_y_minus) * T_y_minus) / two_dy;

            d_pei_dx = kBScaled * ( 2 * get_n_i(u[i+1][j][0], mfi_x_plus) * T_x_plus - 2 * get_n_i(u[i-1][j][0], mfi_x_minus) * T_x_minus) / two_dx;
            d_pei_dy = kBScaled * ( 2 * get_n_i(u[i][j+1][0], mfi_y_plus) * T_y_plus - 2 * get_n_i(u[i][j-1][0], mfi_y_minus) * T_y_minus) / two_dy;

            /* Bx derivatives */
            dBx_dx = ( u[i+1][j][5] - u[i-1][j][5] ) / two_dx;
            dBx_dy = ( u[i][j+1][5] - u[i][j-1][5] ) / two_dy;

            /* Bx derivatives */
            dBy_dx = ( u[i+1][j][6] - u[i-1][j][6] ) / two_dx;
            dBy_dy = ( u[i][j+1][6] - u[i][j-1][6] ) / two_dy;

            /* Bz derivatives */
            dBz_dx = ( u[i+1][j][7] - u[i-1][j][7] ) / two_dx;
            dBz_dy = ( u[i][j+1][7] - u[i][j-1][7] ) / two_dy;

            // COMPOUND TERMS 

            /* j x B terms */
            mag_prefac = mass_frac_n / alpha_in;
            mag_term_x = mag_prefac * (By * dBx_dy - By * dBy_dx - Bz * dBz_dx);
            mag_term_y = mag_prefac * (Bx * dBy_dx - Bx * dBx_dy - Bz * dBz_dy);
            mag_term_z = mag_prefac * (Bx * dBz_dx + By * dBz_dy);

            /* j terms */
            J_term_prefac = alpha_en / (eChargeScaled * n_i * alpha_n);
            J_x = dBz_dy;
            J_y = - dBz_dx;
            J_z = dBy_dx - dBx_dy;

            /* G */
            G_x = mass_frac_n * d_pei_dx - mass_frac_i * d_pn_dx;
            G_y = mass_frac_n * d_pei_dy - mass_frac_i * d_pn_dy;
            G_z = 0;

            /* size comparison */
            // if (i == round(mesh.nCellsX/2) && j == 1){
            //     cout << "1 / a_n: " << 1 / alpha_n << endl;
            //     cout << "G_x term: " << G_x / alpha_n << endl;
            //     cout << "mag term: " << mag_term_x << endl;
            //     cout << "J_z term: " << J_z << endl;
            //     cout << "J_z term * prefac: " << J_term_prefac * J_z << endl; 
            //     cout << "prefac: " << J_term_prefac << endl;
            // }

            // setting and assembling
            u[i][j][9] = - G_x / alpha_n + mag_term_x + J_term_prefac * J_x;
            u[i][j][10] = - G_y / alpha_n + mag_term_y + J_term_prefac * J_y;
            u[i][j][11] = - G_z / alpha_n + mag_term_z + J_term_prefac * J_z;

        }
    }




}






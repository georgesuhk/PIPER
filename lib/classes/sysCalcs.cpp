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


double SysCalcs::get_E(CellVec& u){

    return u[4] - get_MagE(u);  
}


double SysCalcs::get_e(CellVec& u){

    return (get_E(u) - get_KE(u)) / u[0];
};

double SysCalcs::get_p(CellVec& u){
    double e = get_e(u);
    return EoSPtr->interp_p(u[0], e);
};

double SysCalcs::get_p(int i, int j){
    return EoSPtr->get_p(i, j);
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

/* returns the speed of the fastest wave in the system (hyperbolic wave speed) */
double SysCalcs::getFastestWaveSpeed(CellVec& u){

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







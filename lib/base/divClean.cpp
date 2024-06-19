#include "divClean.hpp"

/* Getting parabolic wave speed based on an empirically determined factor */
double getCp(double ch){
    double scalingFactor = 0.18; // empirically determined 
    return sqrt(scalingFactor * ch);
}

/* The parabolic source term update for psi */
double paraUpdatePsi(double psiOld, double ch, double dt){

    double cp = getCp(ch);
    double psiNew = psiOld * exp(-1 * dt * pow(ch, 2.0) / pow(cp, 2.0));
    return psiNew;
}

/* return mid state for Riemann problem encountered in divergence cleaning */
array<double,2> getDCMidState(const CellVec& uL, const CellVec& uR, double ch, char axis){
    double BNorm_L, BNorm_R;
    double psi_L = uL[8];
    double psi_R = uR[8];

    if (axis == 'x'){
        BNorm_L = uL[5];
        BNorm_R = uR[5];
    } else if (axis == 'y'){
        BNorm_L = uL[6];
        BNorm_R = uR[6];
    }

    //the middle value of B norm and psi
    double BNorm_mid = 0.5 * (BNorm_L + BNorm_R) - ( 1 / (2 * ch)) * (psi_R - psi_L);
    double psi_mid = 0.5 * (psi_L + psi_R) - (ch / 2) * (BNorm_R - BNorm_L); 

    return {BNorm_mid, psi_mid};
}

/* return the flux for the separated GLM divergence cleaning system */
array<double,2> getDCFlux(const CellVec& uL, const CellVec& uR, double ch, char axis){

    array<double,2> DCMidState = getDCMidState(uL, uR, ch, axis);
    double BNorm_mid = DCMidState[0];
    double psi_mid = DCMidState[1];

    return {psi_mid, pow(ch, 2.0) * BNorm_mid };
}
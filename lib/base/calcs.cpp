#include "calcs.hpp"

// GENERAL CALCS ======

/* returns temperature according to ideal gas equation */
double getIdealGasTemp(double rho, double p, double particleMass){
    return (particleMass * p) / (rho * kB);
}

/* returns density according to ideal gas equation */
double getIdealGasDensity(double p, double temp, double particleMass){
    return (particleMass * p) / (kB * temp);
}

/* converts temperature from K to eV */
double K_to_eV(double T_K){
    return kB * T_K / elementaryCharge;
}

/* returns the standard reduced mass */
double get_reduced_mass(double m1, double m2){
    return (m1 * m2) / (m1 + m2);
}

/* returns the coulomb logarithm */
double getCoulombLog(double& n_e, double T){
    double output;

    if (T < 580226){ // T < 50eV
        output = 23.4 - 1.15 * log10(n_e / 1e6) + 3.45 * log10(T);
    } else {         // T > 50eV
        output = 25.3 - 1.15 * log10(n_e / 1e6) + 2.3 * log10(T);
    }

    return output;
}

// Collision rates (nu), scaled by sqrt(pAtmos) due to seconds scaling ------

/* returns nu_in, the collision rate between ions and neutrals */
double get_coll_freq_in(double& n_n, double& m_i, double& m_n, double& T){
    double reduced_m_i_n = get_reduced_mass(m_i, m_n);
    double coll_freq = n_n * sqrt( (8 * kB * T) / (myPI * reduced_m_i_n) ) * sigma_coll_in;
    return coll_freq * sqrt(pAtmos);
}

/* returns nu_en, the collision rate between e- and neutrals */
double get_coll_freq_en(double n_n, double& m_i, double& m_n, double& T){
    double reduced_m_i_n = get_reduced_mass(m_i, m_n);
    double coll_freq = n_n * sqrt( (8 * kB * T) / (myPI * reduced_m_i_n) ) * sigma_coll_en;
    return coll_freq * sqrt(pAtmos);
}

/* returns nu_ei, the collision rate between e- and ions */
double get_coll_freq_ei(double& n_i, double& coulombLog, double& T){
    double coll_freq = 3.7e-6 * ( (n_i * coulombLog) / (pow(T, 3.0/2.0)) ) * brag_a0;
    return coll_freq * sqrt(pAtmos);
}



// Collision coefficients (alpha) ------

/* returns the i-n alpha collision coefficients in Braginskii 1965 */
double get_coll_coeff_in(double& n_i, double& n_n, double& m_i, double m_n, double& T){
    double reduced_m_i_n = get_reduced_mass(m_i, m_n);
    return m_i * n_i * get_coll_freq_in(n_n, m_i, m_n, T);
}

/* returns the e-n alpha collision coefficients in Braginskii 1965 */
double get_coll_coeff_en(double& n_i, double& n_n, double& m_i, double& m_n, double& T){
    double reduced_m_i_n = get_reduced_mass(m_i, m_n);
    return m_i * n_i * get_coll_freq_en(n_n, m_i, m_n, T);
}

/* returns the e-i alpha collision coefficients in Braginskii 1965 */
double get_coll_coeff_ei(double n_e, double& n_i, double& coulombLog, double& T){
    return electronMass * n_e * get_coll_freq_ei(n_i, coulombLog, T) * brag_a0;
}



// Collision times (time between collisions) ------ WORK IN PROGRESS

/* returns electron collision time */
// double get_coll_time_e(double Z, double& n_i, double& coulombLog, double& T){
//     double coeff = 2.75e5; // middle ground between Brag and Khomenko
//     // Brag uses 2.8e5, Khomenko uses 2.7e5
//     return coeff * (pow(T, 1.5))/(Z * n_i * coulombLog);
// }




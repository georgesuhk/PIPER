#include "EoSUtils.hpp"


/** 
 * Version of Bisect root finder used when constructing EoS, to find T from rho and p
 * Looks for the value of *T* that will return a given value of *rho* for a specific *p*
*/
double BisectSolverEoS(double& rho, double& p, double& particleMass, Mixture& mix, double atol, double rtol, int maxSteps, double searchTol){

    /* initial guesses for T */
    // current defintion will only work for hydrogen, based on how its average mass per particle changes
    double T_lower = getIdealGasTemp(rho, p, particleMass/(2+searchTol));
    double T_higher = getIdealGasTemp(rho, p, particleMass*(1+searchTol));

    /* midpoint values for T, and the rho obtained using it */
    double T_mid, rho_mid;

    /* difference between rho_mid and target rho */
    double f_mid;

    int step = 0;
    while (fabs(T_lower - T_higher) > atol && fabs(T_lower - T_higher)/T_lower > rtol && step < maxSteps){
        step += 1;

        // compute midpoint
        T_mid = (T_lower + T_higher)/2;

        // find rho value at midpoint
        mix.equilibrate(T_mid, p);
        rho_mid = mix.density();

        // find difference
        f_mid = rho_mid - rho;

        if (f_mid > 0){
            T_lower = T_mid;
        } else {
            T_higher = T_mid;
        }
    }

    return (T_lower + T_higher) / 2;
}

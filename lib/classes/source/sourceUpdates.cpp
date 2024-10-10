#include "sourceUpdates.hpp"

// IMPLICIT ======

/**
 * Ohmic diffusion source term + Joule heating
 * Probs should split into 2 separate functions or at least make joule heating second order
 */
void ohmic_diffusion(Vec2D& u, shared_ptr<SysCalcs> sysPtr, Mesh2D& mesh, double& dt, string diffBCType, BCFunc BC){
    CN_Diffusion_Solver(u, {5,6}, get_Resis, sysPtr, mesh, dt, diffBCType, BC);
}




/**
 * Thermal conduction of the bulk fluid
 */
void conduction(Vec2D& u, shared_ptr<SysCalcs> sysPtr, Mesh2D& mesh, double& dt, string diffBCType, BCFunc BC){
    CN_Conduction_Solver_ADI(u, {4}, get_therm_con, sysPtr, mesh, dt, diffBCType, BC);
}

#include "sourceUpdates.hpp"

// IMPLICIT ======

/**
 * Ohmic diffusion source term + Joule heating
 * Probs should split into 2 separate functions or at least make joule heating second order
 */
void ohmic_diffusion(Vec2D& u, shared_ptr<SysCalcs> sysPtr, Mesh2D& mesh, double& dt, string diffBCType, BCFunc BC){
    CN_Diffusion_Solver(u, 5, get_Resis, sysPtr, mesh, dt, diffBCType, BC);
    CN_Diffusion_Solver(u, 6, get_Resis, sysPtr, mesh, dt, diffBCType, BC);

    //updating total energy
    double resistivity, laplacianBx, laplacianBy;
    double dBz_dx, dBz_dy, dBy_dx, dBx_dy;
    double S;

    for (int i = 1; i < mesh.nCellsX+1; i++){
        for (int j = 1; j < mesh.nCellsY+1; j++){
            resistivity = get_Resis(u[i][j], i, j, sysPtr, true);

            //induction eqn source terms
            laplacianBx = ( u[i+1][j][5] - 2*u[i][j][5] + u[i-1][j][5] )/(mesh.dx) + ( u[i][j+1][5] - 2*u[i][j][5] + u[i][j-1][5] )/(mesh.dy);
            laplacianBy = ( u[i+1][j][6] - 2*u[i][j][6] + u[i-1][j][6] )/(mesh.dx) + ( u[i][j+1][6] - 2*u[i][j][6] + u[i][j-1][6] )/(mesh.dy);

            //energy source terms
            dBz_dx = ( u[i+1][j][7] - u[i-1][j][7] ) / (2*mesh.dx);
            dBz_dy = ( u[i][j+1][7] - u[i][j-1][7] ) / (2*mesh.dy);
            dBy_dx = ( u[i+1][j][6] - u[i-1][j][6] ) / (2*mesh.dx);
            dBx_dy = ( u[i][j+1][5] - u[i][j-1][5] ) / (2*mesh.dy);

            //assembling
            S = resistivity * ( (laplacianBx*u[i][j][5] + laplacianBy*u[i][j][6]) + ( pow(dBz_dy, 2.0) + pow(dBz_dx, 2.0) + pow((dBy_dx - dBx_dy), 2.0) ) );
            u[i][j][4] = u[i][j][4] + S * dt;

        }
    }

    BC(u, mesh, sysPtr);
    sysPtr->getEoSPtr()->cacheAll(u, mesh);
}




/**
 * Thermal conduction of the bulk fluid
 */
void conduction(Vec2D& u, shared_ptr<SysCalcs> sysPtr, Mesh2D& mesh, double& dt, string diffBCType, BCFunc BC){
    CN_Conduction_Solver(u, 4, get_therm_con, sysPtr, mesh, dt, diffBCType, BC);
}

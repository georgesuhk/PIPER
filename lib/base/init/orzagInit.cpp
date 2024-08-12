#include "orzagInit.hpp"


Vec2D initOrzag(vector<double>& initParams, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr, BCFunc BCFunc){
    Vec2D uInit = makeVec2D(mesh.nCellsX + 2, mesh.nCellsY + 2);
    double x, y;
    double gamma = initParams[0];

    for (int i = 1; i < mesh.nCellsX+1; i++){
        for (int j = 1; j < mesh.nCellsY+1; j++){

            x = indexToPos(mesh, i, 'x');
            y = indexToPos(mesh, j, 'y');

            double rho = pow(gamma, 2.0)*rho_SF;
            double vx = -1*sin(2*myPI*y)*sqrt(p_SF)/sqrt(rho_SF);
            double vy = sin(2*myPI*x)*sqrt(p_SF)/sqrt(rho_SF);
            double vz = 0.0*sqrt(p_SF)/sqrt(rho_SF);
            double p = gamma*p_SF;
            double Bx = -1*sin(2*myPI*y)*sqrt(p_SF);
            double By = sin(4*myPI*x)*sqrt(p_SF);
            double Bz = 0.0*sqrt(p_SF);
            double psi = 0.0;
            double wx = 0.0;
            double wy = 0.0;
            double wz = 0.0;

            CellVec uPrim = {rho, vx, vy, vz, p, Bx, By, Bz, psi, wx, wy, wz};

            uInit[i][j] = sysPtr->primToConserv(uPrim);

        }
            
    }
    
    BCFunc(uInit, mesh, sysPtr);
    return uInit;
}

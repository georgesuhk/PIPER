#include "detachInit.hpp"


Vec1D uLeftInit;

Vec2D initDetachment(vector<double>& initParams, string& mixName, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr, BCFunc BCFunc){
    Vec2D uInit = makeVec2D(mesh.nCellsX + 2, mesh.nCellsY + 2);

    // unpacking init params ------
    double p = initParams[0];
    double pSI = p * pAtmos;

    /* upstream temperature */
    double temp_up = initParams[1];

    /* target temperature */
    double temp_target = initParams[2];
    double Bx = initParams[3];
    double By = initParams[4];

    vector<double> temp_range = linspace(temp_up, temp_target, mesh.nCellsX);

    // scaling
    double CsScale = 1.0/sqrt(pAtmos);

    // obtaining bohm boundary conditions
    Mixture mix(mixName);

    double x, y, temp;
    double rho, vx, vy, vz, Bz, psi, wx, wy, wz;
    for (int i = 1; i < mesh.nCellsX+1; i++){
        temp = temp_range[i-1];
        mix.equilibrate(temp, pSI);

        for (int j = 1; j < mesh.nCellsY+1; j++){

            rho = mix.density();
            vx = 0;
            vy = 0;
            vz = 0;
            Bz = 0.0;
            psi = 0.0;
            wx = 0.0;
            wy = 0.0;
            wz = 0.0;

            CellVec uPrim = {rho, vx, vy, vz, p, Bx, By, Bz, psi, wx, wy, wz};

            uInit[i][j] = sysPtr->primToConserv(uPrim);
        }
    }
    uInit[0] = uInit[1];
    uLeftInit = uInit[0];
    BCFunc(uInit, mesh, sysPtr);
    return uInit;
}



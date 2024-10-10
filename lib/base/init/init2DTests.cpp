#include "init2DTests.hpp"


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


Vec2D init_current_sheet(vector<double>& initParams, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr, BCFunc BCFunc){
    Vec2D uInit = makeVec2D(mesh.nCellsX + 2, mesh.nCellsY + 2);
    double x, y;
    double By;
    double v_SF = sqrt(p_SF)/sqrt(rho_SF);

    for (int i = 1; i < mesh.nCellsX+1; i++){
        for (int j = 1; j < mesh.nCellsY+1; j++){
        
            x = indexToPos(mesh, i, 'x');
            y = indexToPos(mesh, j, 'y');

            double rho = 1.0*rho_SF;
            double vx  = 0.1*sin(2*myPI*y)*v_SF;
            double vy  = 0.0;
            double vz  = 0.0;
            double p   = 0.3*p_SF;
            double Bx  = 0.0;
            if (fabs(x) < 0.25){
                By = 1;
            } else {
                By = -1;
            }
            double Bz  = 0.0;
            double psi = 0.0;
            double wx  = 0.0;
            double wy  = 0.0;
            double wz  = 0.0;

            CellVec uPrim = {rho, vx, vy, vz, p, Bx, By, Bz, psi, wx, wy, wz};

            uInit[i][j] = sysPtr->primToConserv(uPrim);

        }
            
    }
    
    BCFunc(uInit, mesh, sysPtr);
    return uInit;
}


Vec2D init_reconnection(vector<double>& initParams, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr, BCFunc BCFunc){
    Vec2D uInit = makeVec2D(mesh.nCellsX + 2, mesh.nCellsY + 2);
    double x, y;
    double By, Bz;
    double rho, p;
    double Lr = 0.05;

    for (int i = 1; i < mesh.nCellsX+1; i++){
        for (int j = 1; j < mesh.nCellsY+1; j++){
        
            x = indexToPos(mesh, i, 'x');
            y = indexToPos(mesh, j, 'y');

            p = 1.0 * p_SF;
            if (fabs(x) < Lr && fabs(y) < 2 * Lr){
                // rho = 15 * rho_SF * (cos(myPI * x / (2*Lr))+1);
                rho = 8 * rho_SF * (cos(myPI * x / (2*Lr))+1);

            } else {
                rho = 1.0 * rho_SF;
            }
            double vx  = 0.0;
            double vy  = 0.0;
            double vz  = 0.0;
            double Bx  = 0.0;
            if (x < -Lr){
                By = -1*sqrt(p_SF);
                Bz = 0*sqrt(p_SF);
            } else if (x <= Lr) { 
                By = sin(myPI * x * Lr / 2)*sqrt(p_SF);
                Bz = cos(myPI * x * Lr / 2)*sqrt(p_SF);
            } else {
                By = 1*sqrt(p_SF);
                Bz = 0*sqrt(p_SF);
            }
            double psi = 0.0;
            double wx  = 0.0;
            double wy  = 0.0;
            double wz  = 0.0;

            CellVec uPrim = {rho, vx, vy, vz, p, Bx, By, Bz, psi, wx, wy, wz};

            uInit[i][j] = sysPtr->primToConserv(uPrim);

        }
            
    }
    
    BCFunc(uInit, mesh, sysPtr);
    return uInit;
}

Vec2D init_mag_recon2(vector<double>& initParams, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr, BCFunc BCFunc){
    Vec2D uInit = makeVec2D(mesh.nCellsX + 2, mesh.nCellsY + 2);
    double x, y;
    double By, Bz;
    double rho, p;
    double Lr = 0.05;

    for (int i = 1; i < mesh.nCellsX+1; i++){
        for (int j = 1; j < mesh.nCellsY+1; j++){
        
            x = indexToPos(mesh, i, 'x');
            y = indexToPos(mesh, j, 'y');

            p = 1.0 * p_SF;
            if (fabs(y) < 2 * Lr){
                if (fabs(x) < Lr){
                    if (x < 0){
                        rho = 1.6 * rho_SF;
                    } else {
                        rho = 1.6 * rho_SF;
                    }
                } else {
                    rho = 1.0 * rho_SF;
                }
            } else {
                rho = 1.0 * rho_SF;
            }

            double vx  = 0.0;
            double vy  = 0.0;
            double vz  = 0.0;
            double Bx  = 0.0;
            if (x < -Lr){
                By = -1*sqrt(p_SF);
                Bz = 0*sqrt(p_SF);
            } else if (x <= Lr) { 
                By = sin(myPI * x * Lr / 2)*sqrt(p_SF);
                Bz = cos(myPI * x * Lr / 2)*sqrt(p_SF);
            } else {
                By = 1*sqrt(p_SF);
                Bz = 0*sqrt(p_SF);
            }
            double psi = 0.0;
            double wx  = 0.0;
            double wy  = 0.0;
            double wz  = 0.0;

            CellVec uPrim = {rho, vx, vy, vz, p, Bx, By, Bz, psi, wx, wy, wz};

            uInit[i][j] = sysPtr->primToConserv(uPrim);

        }
            
    }
    
    BCFunc(uInit, mesh, sysPtr);
    return uInit;
}





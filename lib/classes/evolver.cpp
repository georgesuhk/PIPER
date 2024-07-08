#include "evolver.hpp"

// EVOLVER BASE CLASS FUNCTIONS ======

void Evolver::updateCh(Vec2D& u, shared_ptr<SysCalcs> sysPtr, Mesh2D& mesh){
    double cMax = 0;
    double c;

    //Looping over all cells to find maximum
    for (int i = 1; i < mesh.nCellsX+1; i++){
        for (int j = 1; j < mesh.nCellsY+1; j++){

            //get fastest wave speed for cell
            c = sysPtr->getFastestWaveSpeed(u[i][j]);
            if (c > cMax){ cMax = c;}
        }
    }

    ch = cMax;
}

double Evolver::getTimeStep(Mesh2D& mesh){
    double minCellGap = minBetween(mesh.dx, mesh.dy);
    // cout << "ch: " << ch << endl;
    return Cnum * (minCellGap / ch);
}

void Evolver::evolveMat(Vec2D& u, shared_ptr<SysCalcs> sysPtr, Mesh2D& mesh, double& dt, char axis){
    Vec2D uNext = makeVec2D(mesh.nCellsX+2, mesh.nCellsY+2);

    //x update ------------
    if (axis == 'x'){
        updateFlux(u, sysPtr, mesh, dt, 'x');

        for (int i = 1; i < mesh.nCellsX+1; i++){
            for (int j = 1; j < mesh.nCellsY+1; j++){

                // updating each variable
                for (int var = 0; var < (cellVarsNums); var++){
                    uNext[i][j][var] = u[i][j][var] - (dt/mesh.dx) * (flux[i+1][j][var] - flux[i][j][var]);
                }

                // parabolic divergence cleaning term
                if (doDC){
                    uNext[i][j][8] = paraUpdatePsi(uNext[i][j][8], ch, dt);
                }
            }
        }

    //Y update -------------
    } else if (axis == 'y'){
        updateFlux(u, sysPtr, mesh, dt, 'y');

        for (int i = 1; i < mesh.nCellsX+1; i++){
            for (int j = 1; j < mesh.nCellsY+1; j++){

                // updating each variable
                for (int var = 0; var < (cellVarsNums); var++){
                    uNext[i][j][var] = u[i][j][var] - (dt/mesh.dy) * (flux[i][j+1][var] - flux[i][j][var]);
                }

                // parabolic divergence cleaning term
                if (doDC){
                    uNext[i][j][8] = paraUpdatePsi(uNext[i][j][8], ch, dt);
                }
            }
        }
    }

    evolverBCFunc(uNext, mesh, sysPtr);
    u = uNext;
}

Vec2D Evolver::getFlux(Vec2D& u, shared_ptr<SysCalcs> sysPtr, Mesh2D& mesh, double& dt, char axis){
    // check if flux object has been initialised
    if (flux.empty()){
        updateFlux(u, sysPtr, mesh, dt, axis);        
    }

    return flux;
}

void Evolver::updateFlux(Vec2D& u, shared_ptr<SysCalcs> sysPtr, Mesh2D& mesh, double& dt, char axis){
    throw runtime_error("updateFlux() should not be called in base Evolver class");
}

string Evolver::getEvolverType(){
    return evolverType;
}

void Evolver::setDoDC(bool inputDoDC){
    doDC = inputDoDC;
}

bool Evolver::getDoDC(){
    return doDC;
}



// FORCE Evolver functions ======

void FORCEEvolver::updateFlux(Vec2D& u, shared_ptr<SysCalcs> sysPtr, Mesh2D& mesh, double& dt, char axis) {
    Vec2D fluxFORCE = makeVec2D(mesh.nCellsX+2, mesh.nCellsY+2);

    //calculating component (LF and RI) fluxes
    Vec2D fluxLF = getFluxLF(u, sysPtr, mesh, dt, axis);
    Vec2D fluxRI = getFluxRI(u, sysPtr, mesh, dt, axis);
    array<double, 2> DCFlux;

    //combining to make FORCE
    for (int i = 1; i < mesh.nCellsX+2; i++){
        for (int j = 1; j < mesh.nCellsY+2; j++){
            for (int var = 0; var < (cellVarsNums); var++){
                fluxFORCE[i][j][var] = 0.5 * ( fluxLF[i][j][var] + fluxRI[i][j][var] );

                // divergence cleaning functionality
                if (doDC){
                    if (axis == 'x'){
                        DCFlux = getDCFlux(u[i-1][j], u[i][j], ch, axis);
                        fluxFORCE[i][j][5] = DCFlux[0];
                        fluxFORCE[i][j][8] = DCFlux[1];          
                    } else if (axis == 'y'){
                        DCFlux = getDCFlux(u[i][j-1], u[i][j], ch, axis);
                        fluxFORCE[i][j][6] = DCFlux[0];
                        fluxFORCE[i][j][8] = DCFlux[1];  
                    }
                }
            }
        }
    }

    flux = fluxFORCE;
}




// SLIC Evolver Functions ======

void SLICEvolver::setW(double inputW){
    w = inputW;
}

double SLICEvolver::getW(){
    return w;
}

void SLICEvolver::setLimiter(Limiter inputLimFunc){
    limFunc = inputLimFunc;
}

Limiter SLICEvolver::getLimiter(){
    return limFunc;
}

void SLICEvolver::updateFlux(Vec2D& u, shared_ptr<SysCalcs> sysPtr, Mesh2D& mesh, double& dt, char axis){

    //reconstruction
    array<Vec2D,2> uRecon = slopeRecon(u, mesh, w, axis, sysPtr, limFunc, evolverBCFunc);

    //local half step update
    array<Vec2D,2> uHS = reconHS(uRecon[0], uRecon[1], sysPtr, mesh, dt, axis, evolverBCFunc);

    //calculating component (LF and RI) fluxes
    Vec2D fluxLF = getFluxLF(uHS[0], uHS[1], sysPtr, mesh, dt, axis);
    Vec2D fluxRI = getFluxRI(uHS[0], uHS[1], sysPtr, mesh, dt, axis);

    //combining to make SLIC
    Vec2D fluxSLIC = makeVec2D(mesh.nCellsX+2, mesh.nCellsY+2);

    array<double, 2> DCFlux;
    for (int i = 1; i < mesh.nCellsX+2; i++){
        for (int j = 1; j < mesh.nCellsY+2; j++){
            for (int var = 0; var < (cellVarsNums); var++){
                fluxSLIC[i][j][var] = 0.5 * ( fluxLF[i][j][var] + fluxRI[i][j][var] );

                // divergence cleaning functionality
                if (doDC){
                    if (axis == 'x'){
                        DCFlux = getDCFlux(u[i-1][j], u[i][j], ch, axis);
                        fluxSLIC[i][j][5] = DCFlux[0];
                        fluxSLIC[i][j][8] = DCFlux[1];          
                    } else if (axis == 'y'){
                        DCFlux = getDCFlux(u[i][j-1], u[i][j], ch, axis);
                        fluxSLIC[i][j][6] = DCFlux[0];
                        fluxSLIC[i][j][8] = DCFlux[1];  
                    }
                }
            }
        }
    }

    flux = fluxSLIC;
}
#include "fluxCalcs.hpp"

/* LF Flux*/
Vec2D getFluxLF(Vec2D& u, shared_ptr<SysCalcs> sysPtr, const Mesh2D& mesh, double dt, char axis){
    Vec2D flux = makeVec2D(mesh.nCellsX+2, mesh.nCellsY+2);

    //x direction
    if (axis == 'x'){
        for (int i = 1; i < mesh.nCellsX+2; i++){
            for (int j = 1; j < mesh.nCellsY+2; j++){
                for (int var = 0; var < (cellVarsNums); var++){
                    flux[i][j][var] = 0.5 * (mesh.dx/dt) * (u[i-1][j][var] - u[i][j][var]) + 0.5 * ( sysPtr->f(u[i][j], i, j)[var] + sysPtr->f(u[i-1][j], i-1, j)[var] );
                }
            }
        }
    }

    //y direction
    else if (axis == 'y'){
        for (int i = 1; i < mesh.nCellsX+2; i++){
            for (int j = 1; j < mesh.nCellsY+2; j++){
                for (int var = 0; var < (cellVarsNums); var++){
                    flux[i][j][var] = 0.5 * (mesh.dy/dt) * (u[i][j-1][var] - u[i][j][var]) + 0.5 * ( sysPtr->g(u[i][j], i, j)[var] + sysPtr->g(u[i][j-1], i, j-1)[var]);
                }
            }
        }
    }

    else {
        throw invalid_argument("In getfluxLF(): axis given is invalid");
    }

    return flux;
}

/* LF for SLIC */
Vec2D getFluxLF(Vec2D& uLHS, Vec2D& uRHS, shared_ptr<SysCalcs> sysPtr, const Mesh2D& mesh, double dt, char axis){
    Vec2D flux = makeVec2D(mesh.nCellsX+2, mesh.nCellsY+2);

    //x direction
    if (axis == 'x'){
        for (int i = 1; i < mesh.nCellsX+2; i++){
            for (int j = 1; j < mesh.nCellsY+2; j++){
                for (int var = 0; var < (cellVarsNums); var++){
                    flux[i][j][var] = 0.5 * (mesh.dx/dt) * (uRHS[i-1][j][var] - uLHS[i][j][var]) + 0.5 * ( sysPtr->f(uLHS[i][j], i, j, true)[var] + sysPtr->f(uRHS[i-1][j], i-1, j, true)[var] );
                }
            }
        }
    }

    //y direction
    else if (axis == 'y'){
        for (int i = 1; i < mesh.nCellsX+2; i++){
            for (int j = 1; j < mesh.nCellsY+2; j++){
                for (int var = 0; var < (cellVarsNums); var++){
                    flux[i][j][var] = 0.5 * (mesh.dy/dt) * (uRHS[i][j-1][var] - uLHS[i][j][var]) + 0.5 * ( sysPtr->g(uLHS[i][j], i, j, true)[var] + sysPtr->g(uRHS[i][j-1], i, j-1, true)[var]);
                }
            }
        }
    }

    else {
        throw invalid_argument("In getfluxLF(): axis given is invalid");
    }

    return flux;
}

/* RI Flux */
Vec2D getFluxRI(Vec2D& u, shared_ptr<SysCalcs> sysPtr, const Mesh2D& mesh, double dt, char axis){
    Vec2D flux = makeVec2D(mesh.nCellsX+2, mesh.nCellsY+2);
    Vec2D uHalfStepped = makeVec2D(mesh.nCellsX+2, mesh.nCellsY+2);

    //x direction
    if (axis == 'x'){
        for (int i = 1; i < mesh.nCellsX+2; i++){
            for (int j = 1; j < mesh.nCellsY+2; j++){
                for (int var = 0; var < (cellVarsNums); var++){
                    uHalfStepped[i][j][var] = 0.5 * (u[i-1][j][var] + u[i][j][var]) - 0.5 * (dt/mesh.dx) * ( sysPtr->f(u[i][j], i, j)[var] - sysPtr->f(u[i-1][j], i-1, j)[var] );
                }
            }
        }

        for (int i = 1; i < mesh.nCellsX+2; i++){
            for (int j = 1; j < mesh.nCellsY+2; j++){
                flux[i][j] = sysPtr->f(uHalfStepped[i][j], i, j, true);
            }
        }
    }
    //y direction
    else if (axis == 'y'){
        for (int i = 1; i < mesh.nCellsX+2; i++){
            for (int j = 1; j < mesh.nCellsY+2; j++){
                for (int var = 0; var < (cellVarsNums); var++){
                    uHalfStepped[i][j][var] = 0.5 * (u[i][j-1][var] + u[i][j][var]) - 0.5 * (dt/mesh.dy) * ( sysPtr->g(u[i][j], i, j)[var] - sysPtr->g(u[i][j-1], i, j-1)[var] );
                }
            }
        }

        for (int i = 1; i < mesh.nCellsX+2; i++){
            for (int j = 1; j < mesh.nCellsY+2; j++){
                flux[i][j] = sysPtr->g(uHalfStepped[i][j], i, j, true);
            }
        }
    }
    else {
        throw invalid_argument("In getfluxRI(): axis given is invalid");
    }

    return flux;
}

// RI for SLIC
Vec2D getFluxRI(Vec2D& uLHS, Vec2D& uRHS, shared_ptr<SysCalcs> sysPtr, const Mesh2D& mesh, double dt, char axis){
    Vec2D flux = makeVec2D(mesh.nCellsX+2, mesh.nCellsY+2);
    Vec2D uHalfStepped = makeVec2D(mesh.nCellsX+2, mesh.nCellsY+2);

    //x direction
    if (axis == 'x'){
        for (int i = 1; i < mesh.nCellsX+2; i++){
            for (int j = 1; j < mesh.nCellsY+2; j++){
                for (int var = 0; var < (cellVarsNums); var++){
                    uHalfStepped[i][j][var] = 0.5 * (uRHS[i-1][j][var] + uLHS[i][j][var]) - 0.5 * (dt/mesh.dx) * ( sysPtr->f(uLHS[i][j], i, j, true)[var] - sysPtr->f(uRHS[i-1][j], i-1, j, true)[var] );
                }
            }
        }

        for (int i = 1; i < mesh.nCellsX+2; i++){
            for (int j = 1; j < mesh.nCellsY+2; j++){
                flux[i][j] = sysPtr->f(uHalfStepped[i][j], i, j, true);
            }
        }
    }
    //y direction
    else if (axis == 'y'){
        for (int i = 1; i < mesh.nCellsX+2; i++){
            for (int j = 1; j < mesh.nCellsY+2; j++){
                for (int var = 0; var < (cellVarsNums); var++){
                    uHalfStepped[i][j][var] = 0.5 * (uRHS[i][j-1][var] + uLHS[i][j][var]) - 0.5 * (dt/mesh.dy) * ( sysPtr->g(uLHS[i][j], i, j, true)[var] - sysPtr->g(uRHS[i][j-1], i, j-1, true)[var] );
                }
            }
        }

        for (int i = 1; i < mesh.nCellsX+2; i++){
            for (int j = 1; j < mesh.nCellsY+2; j++){
                flux[i][j] = sysPtr->g(uHalfStepped[i][j], i, j, true);
            }
        }
    }
    else {
        throw invalid_argument("In getfluxRI(): axis given is invalid");
    }

    return flux;
}



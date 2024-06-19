#include "reconstruct.hpp"

Vec2D getDeltas(const Vec2D& u, Mesh2D& mesh, double w, char axis){
    Vec2D deltaVec = makeVec2D(mesh.nCellsX+1, mesh.nCellsY+1); //delta vec will have a ghost cells at the start that should never be accessed
    double deltaMinus, deltaPlus;

    if (axis == 'x'){
        for (int i = 1; i < mesh.nCellsX+1; i++){
            for (int j = 1; j < mesh.nCellsY+1; j++){
                for (int var = 0; var < (cellVarsNums); var++){
                    deltaMinus = u[i][j][var] - u[i-1][j][var];
                    deltaPlus = u[i+1][j][var] - u[i][j][var];

                    deltaVec[i][j][var] = 0.5 * (1 + w) * deltaMinus + 0.5 * (1 - w) * deltaPlus; 
                }

                //no reconstruct for psi
                deltaVec[i][j][8] = 0;
            }
        }
    } else if (axis == 'y'){
        for (int i = 1; i < mesh.nCellsX+1; i++){
            for (int j = 1; j < mesh.nCellsY+1; j++){
                for (int var = 0; var < (cellVarsNums); var++){
                    deltaMinus = u[i][j][var] - u[i][j-1][var];
                    deltaPlus = u[i][j+1][var] - u[i][j][var];

                    deltaVec[i][j][var] = 0.5 * (1 + w) * deltaMinus + 0.5 * (1 - w) * deltaPlus; 
                }

                //no reconstruct for psi
                deltaVec[i][j][8] = 0;
            }
        }
    } else {
        throw invalid_argument("In getDeltas(): axis is invalid");
    }

    return deltaVec;
}


/* Returns slope reconstructed variable vectors uL and uR (in that order) */
array<Vec2D,2> slopeRecon(const Vec2D& u, Mesh2D& mesh, double w, char axis, shared_ptr<SysCalcs> sysPtr, Limiter limFunc, BCFunc BCFunc){
    Vec2D deltaVec = getDeltas(u, mesh, w, axis);
    Vec2D limiterVec = getLimiterVec(u, mesh, limFunc, axis);

    Vec2D uReconL = makeVec2D(mesh.nCellsX+2, mesh.nCellsY+2); //same size as u, need to apply BCs to
    Vec2D uReconR = makeVec2D(mesh.nCellsX+2, mesh.nCellsY+2); //same size as u, need to apply BCs to

    double reconVar;
    for (int i = 1; i < mesh.nCellsX+1; i++){
        for (int j = 1; j < mesh.nCellsY+1; j++){
            for (int var = 0; var < (cellVarsNums); var++){
                reconVar = 0.5 * (limiterVec[i][j][var] * deltaVec[i][j][var]);
                uReconL[i][j][var] = u[i][j][var] - reconVar;
                uReconR[i][j][var] = u[i][j][var] + reconVar;
            }
        }
    }

    BCFunc(uReconL, mesh, sysPtr);
    BCFunc(uReconR, mesh, sysPtr);

    return {uReconL, uReconR};
}


/* Returns half stepped reconstructed states in order of left, right */
array<Vec2D,2> reconHS(Vec2D& uL, Vec2D& uR, shared_ptr<SysCalcs> sysPtr, Mesh2D& mesh, double dt, char axis, BCFunc BCFunc){
    Vec2D uReconHS_L = makeVec2D(mesh.nCellsX+2, mesh.nCellsY+2); //u reconstructed half stepped
    Vec2D uReconHS_R = makeVec2D(mesh.nCellsX+2, mesh.nCellsY+2); //u reconstructed half stepped

    double reconHS_var;
    for (int i = 1; i < mesh.nCellsX+1; i++){
        for (int j = 1; j < mesh.nCellsY+1; j++){

            if (axis == 'x'){
                for (int var = 0; var < (cellVarsNums); var++){
                    reconHS_var = 0.5 * (dt/mesh.dx) * (sysPtr->f(uR[i][j])[var] - sysPtr->f(uL[i][j])[var]);

                    uReconHS_L[i][j][var] = uL[i][j][var] - reconHS_var;
                    uReconHS_R[i][j][var] = uR[i][j][var] - reconHS_var;
                }

            } else if (axis == 'y'){
                for (int var = 0; var < (cellVarsNums); var++){
                    reconHS_var = 0.5 * (dt/mesh.dy) * (sysPtr->g(uR[i][j])[var] - sysPtr->g(uL[i][j])[var]);

                    uReconHS_L[i][j][var] = uL[i][j][var] - reconHS_var;
                    uReconHS_R[i][j][var] = uR[i][j][var] - reconHS_var;
                }
            } else {
                throw invalid_argument("in reconHS(): axis argument invalid");
            }
        }
    }

    BCFunc(uReconHS_L, mesh, sysPtr);
    BCFunc(uReconHS_R, mesh, sysPtr);

    return {uReconHS_L, uReconHS_R};
}
#include "explicitSolvers.hpp"


void RK2(Vec2D& u, shared_ptr<SysCalcs> sysPtr, Mesh2D& mesh, double& dt, SourceFuncEx sourceFunc, BCFunc BC){

    Vec2D k1 = sourceFunc(u, mesh, sysPtr, BC);
    Vec2D inner = makeVec2D(mesh.nCellsX + 2, mesh.nCellsY + 2); //a temp TopVec object for holding data
    
    //making k2 term
    for (int i = 1; i < mesh.nCellsX+1; i++){
        for (int j = 1; j < mesh.nCellsY+1; j++){
            for (int var = 0; var < cellVarsNums; var++){
                inner[i][j][var] = u[i][j][var] + dt * k1[i][j][var];
            }
        }
    }

    BC(inner, mesh, sysPtr);
    Vec2D k2 = sourceFunc(inner, mesh, sysPtr, BC);

    //making uNext
    for (int i = 1; i < mesh.nCellsX+1; i++){
        for (int j = 1; j < mesh.nCellsY+1; j++){
            for (int var = 0; var < cellVarsNums; var++){
                u[i][j][var] = u[i][j][var] + 0.5 * dt * (k1[i][j][var] + k2[i][j][var]);
            }
        }
    }

    BC(u, mesh, sysPtr);
}
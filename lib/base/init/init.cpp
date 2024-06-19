#include "init.hpp"

/**
 * internal function used in init functions that handles the comparison between the current position
 * and the interface intervals specified
*/
int getFitIntNum(double r, vector<double>& interfacePositions){
    int fitIntNum = 0;

    for (int k = 0; k < int(interfacePositions.size()); k++){
        if (r > interfacePositions[k]){
            fitIntNum += 1;
        }     
    }
    return fitIntNum;
}

Vec2D initPlanar(vector<CellVec> initCellVecs, vector<double> interfacePositions, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr, BCFunc BCFunc, char axis){
    Vec2D uInit = makeVec2D(mesh.nCellsX + 2, mesh.nCellsY + 2);

    double r; //position
    for (int i = 1; i < mesh.nCellsX+1; i++){
        for (int j = 1; j < mesh.nCellsY+1; j++){
            if (axis == 'x'){
                r = indexToPos(mesh, i, axis);
            } else if (axis == 'y'){
                r = indexToPos(mesh, j, axis);
            } else {
                throw invalid_argument("invalid axis encountered in initPlanar()");
            }

            //comparing r with interface positions to find which interval it fits into
            int fitIntNum = getFitIntNum(r, interfacePositions);
            uInit[i][j] = sysPtr->primToConserv(initCellVecs[fitIntNum]);
        } 
    }

    // apply boundary conditions
    BCFunc(uInit, mesh, sysPtr);

    // can also cache resistivities here

    return uInit;
}








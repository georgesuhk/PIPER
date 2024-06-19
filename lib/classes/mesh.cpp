#include "mesh.hpp"

double indexToPos(Mesh2D mesh, int idx, char axis){
    double pos;

    if (axis == 'x'){
        pos = mesh.xMin + (idx - 0.5) * mesh.dx;
    } else if (axis == 'y'){
        pos = mesh.yMin + (idx - 0.5) * mesh.dy;
    }

    return pos;

}


int posToIndex(Mesh2D mesh, double pos, char axis){
    int idx;

    if (axis == 'x'){
        idx = round( (pos - mesh.xMin)/mesh.dx + 0.5 );
    } else if (axis == 'y'){
        idx = round( (pos - mesh.yMin)/mesh.dy + 0.5 );
    } else {
        throw invalid_argument("in posToIndex: invalid axis");
    }

    return idx;

}
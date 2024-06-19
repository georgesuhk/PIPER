#ifndef MESH_HPP
#define MESH_HPP


#include "settings.hpp"

struct Mesh2D {
    string meshType = "Mesh2D";

    double xMin;
    double xMax;
    double dx;
    int nCellsX;

    double yMin;
    double yMax;
    double dy;
    int nCellsY;

    //default constructor 
    Mesh2D(){};

    //parameterized (real) constructor
    Mesh2D(double inputXMin, double inputXMax, double inputNCellsX,
    double inputYMin, double inputYMax, double inputNCellsY): 
    xMin(inputXMin), xMax(inputXMax), nCellsX(inputNCellsX), 
    yMin(inputYMin), yMax(inputYMax), nCellsY(inputNCellsY){
        dx = (xMax - xMin)/nCellsX;
        dy = (yMax - yMin)/nCellsY; 
    }
};

double indexToPos(Mesh2D mesh, int idx, char axis);
int posToIndex(Mesh2D mesh, double pos, char axis);

#endif
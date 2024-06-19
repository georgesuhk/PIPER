#ifndef UTILS_HPP
#define UTILS_HPP

#include "settings.hpp"
#include "mesh.hpp"
#include "opOverload.hpp"
#include "mathsUtils.hpp"


// GENERAL MISC UTILS ======

//functions for creating data storage objects of given length
Vec2D makeVec2D(int lenX, int lenY);
Scalar2D makeScalar2D(int lenX, int lenY);


//functions for accessing and manipulating matrices
void setZeros(Vec2D& u);
Scalar1D getColumn(const Scalar2D& matrix, int idx);
Scalar1D getColumnEoS(const Scalar2D& matrix, int idx);
Scalar1D getRow(const Scalar2D& matrix, int idx);
Vec1D getRow(const Vec2D& matrix, int idx);


/* Creates a equally spaced list, similar to range() in python */
vector<int> range(int start, int end, int gap);

/* Similar to python np.linspace() */
vector<double> linspace(double start, double end, int nPoints);

/* generates a table (vector<vector<double>>) from csv data */
vector<vector<double>> table_from_csv(string filename, char delimiter);



#endif
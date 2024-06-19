#ifndef TYPEDEFS_HPP
#define TYPEDEFS_HPP

#include "modules.hpp"

// storage
typedef vector<double> CellVec;  // container for variables within each cell
typedef vector<CellVec> Vec1D;   // container for cellVecs along 1 dimension
typedef vector<Vec1D> Vec2D;     // container for cellVecs in 2 dimensions, a vector of vec1Ds

typedef vector<double> Scalar1D; // container for scalar objects along 1 dimension
typedef vector<Scalar1D> Scalar2D; // container for scalar objects in 2 dimensions, a vector of scalar1Ds

// others
typedef array<int,2> IdxPair; //a pair of indices, usually for index of the cell left and right of an interface

#endif 
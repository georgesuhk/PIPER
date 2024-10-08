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

/* Similar to python np.linspace() but samples expoentially */
vector<double> linspaceLog(double start, double end, int nPoints);

/* calculates current at a given i j position */
Eigen::Vector3d get_J(Vec2D& u, Mesh2D& mesh, int i, int j);

void enforce_symmetry(Vec2D& u, Mesh2D& mesh, double r_tol_x, double r_tol_y, bool enforce_x, bool enforce_y);



#endif
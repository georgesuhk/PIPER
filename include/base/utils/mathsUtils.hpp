#ifndef MATHSUTILS_HPP
#define MATHSUTILS_HPP

#include "settings.hpp"
#include "mesh.hpp"

extern const int cellVarsNums;
// LINEAR ALGEBRA TOOLS ======

array<double, 2> vectorAdd(array<double,2> a, array<double,2> b);
array<double, 3> vectorAdd(array<double,3> a, array<double,3> b);
array<double, 2> vectorSubtract(array<double,2> a, array<double,2> b);
array<double, 3> vectorSubtract(array<double,3> a, array<double,3> b);
array<double,2> scalarTimesVec(double scalar, array<double,2> vec);
array<double,2> vecDivideByScalar(double scalar, array<double,2> vec);
array<double,2> normaliseVec(array<double,2> vec);

/* returns the dot product between 2 vectors */
double dotProduct(array<double, 3> a, array<double, 3> b);
double dotProduct(array<double, 2> a, array<double, 2> b);

/* Returns the magnitude of a vector */
double getVectorMagnitude(array<double, 3> a);
double getVectorMagnitude(array<double, 2> a);

/**
 * Returns the tangent vector to a normal
 * (direction should be assumed to be random)
*/
array<double,2> getTangVec(array<double, 2> normVec);

/* returns closets distance of point to a line */
double distanceToLine(double x, double y, double diagOffset);

/* Thomas algorithm, used to solve tridiagonal matrices*/
vector<double> thomasAlgo(Mesh2D& mesh, vector<double> A_diag, vector<double> A_offdiag_l, vector<double> A_offdiag_r, vector<double> Bvector, char axis);




// OTHER MATHS TOOLS ======

// Simple functions ------

/* Gets the sign of an input value */
int sgn(double x);
int sgn(int x);

double minBetween(double a, double b);
int minBetween(int a, int b);

double maxBetween(double a, double b);
int maxBetween(int a, int b);

/* checks to see if 2 input float values are within a given tolerance */
bool floatsAreClose(double a, double b, double tol);



// Advanced functions ------
/* bilinear interpolation function */
CellVec bilinearInterp(const Vec2D& u, Mesh2D& mesh, double xPos, double yPos);

/* bilinear interpolation function for scalar fields (used in tabulated EoS) */
double bilinearInterp(vector<vector<double>>& data, double& xVar, double& yVar, double& x1, double& x2, double xLowerIdx,
double& y1, double& y2, double yLowerIdx);

/* Uses Bisection search to return lower index that bounds a value in an array, used in TabEoS */
int getLowerBound(const double& val, array<int,2> activeRange, vector<double>& data);

/* Bisection root finder used in TabEoS to find p from rho and e */
double BisectSolver(double& rho, int& rhoLowerIdx, double& e, Scalar1D& rhoData, Scalar2D& eData, Scalar1D pData, double atol, int maxSteps);


#endif
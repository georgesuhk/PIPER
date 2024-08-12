#ifndef IMPLICITSOLVERS_HPP
#define IMPLICITSOLVERS_HPP

#include "utils.hpp"
#include "BC.hpp"

// DIFFUSION COEFFICIENTS ======

/* type def for general diffusion coefficient returning functions */

typedef double (*GetCoeffFunc)(CellVec&, int, int, shared_ptr<SysCalcs>, bool);

/* Returns the resistivity */
double get_Resis(CellVec& u, int i, int j, shared_ptr<SysCalcs> sysPtr, bool interp);



// IMPLICIT SOLVERS ======

/* Crank-Nicolson implicit solver for the diffusion equation */
void CN_Diffusion_Solver(Vec2D& u, int varIdx, GetCoeffFunc getCoeffFunc, shared_ptr<SysCalcs> sysPtr, Mesh2D& mesh, double& dt, string diffBCType, BCFunc BCFunc);

#endif
#ifndef IMPLICITSOLVERS_HPP
#define IMPLICITSOLVERS_HPP

#include "utils.hpp"
#include "BC.hpp"

// DIFFUSION COEFFICIENTS ======

/* type def for general diffusion coefficient returning functions */

typedef double (*GetCoeffFunc)(CellVec&, int, int, shared_ptr<SysCalcs>, bool);

/* Returns the resistivity */
double get_Resis(CellVec& u, int i, int j, shared_ptr<SysCalcs> sysPtr, bool interp);

/* Returns the conductivity */
double get_therm_con(CellVec& u, int i, int j, shared_ptr<SysCalcs> sysPtr, bool interp);


// IMPLICIT SOLVERS ======

/* Non-ADI Crank-Nicolson implicit solver for the diffusion equation */
void CN_Diffusion_Solver(Vec2D& u, vector<int> varIdxList, GetCoeffFunc getCoeffFunc, shared_ptr<SysCalcs> sysPtr, Mesh2D& mesh, double& dt, string diffBCType, BCFunc BCFunc);

/* Crank-Nicolson implicit solver for the diffusion equation */
void CN_Diffusion_Solver_ADI(Vec2D& u, vector<int> varIdxList, GetCoeffFunc getCoeffFunc, shared_ptr<SysCalcs> sysPtr, Mesh2D& mesh, double& dt, string diffBCType, BCFunc BCFunc);

/* Crank-Nicolson implicit solver for thermal conduction */
void CN_Conduction_Solver_ADI(Vec2D& u, vector<int> varIdxList, GetCoeffFunc getCoeffFunc, shared_ptr<SysCalcs> sysPtr, Mesh2D& mesh, double& dt, string diffBCType, BCFunc BCFunc);

#endif
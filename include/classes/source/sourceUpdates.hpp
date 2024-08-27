#ifndef SOURCEUPDATES_HPP
#define SOURCEUPDATES_HPP

#include "implicitSolvers.hpp"
#include "explicitSolvers.hpp"

typedef void (*implicitSource)(Vec2D&, shared_ptr<SysCalcs>, Mesh2D&, double&, string, BCFunc);


void ohmic_diffusion(Vec2D& u, shared_ptr<SysCalcs> sysPtr, Mesh2D& mesh, double& dt, string diffBCType, BCFunc BC);

/* Thermal Conduction */
void conduction(Vec2D& u, shared_ptr<SysCalcs> sysPtr, Mesh2D& mesh, double& dt, string diffBCType, BCFunc BC);



#endif
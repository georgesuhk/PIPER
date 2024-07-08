#ifndef EXPLICITSOLVERS_HPP
#define EXPLICITSOLVERS_HPP

#include "sourceFuncEx.hpp"
typedef void (*ExplicitSolver)(Vec2D&, shared_ptr<SysCalcs>, Mesh2D&, double&, SourceFuncEx, BCFunc);


void RK2(Vec2D& u, shared_ptr<SysCalcs> sysPtr, Mesh2D& mesh, double& dt, SourceFuncEx sourceFunc, BCFunc BC);



#endif 
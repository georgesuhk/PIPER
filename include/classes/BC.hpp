#ifndef BC_HPP
#define BC_HPP

#include "sysCalcs.hpp"
#include "mesh.hpp"

typedef void (*BCFunc)(Vec2D&, Mesh2D&, shared_ptr<SysCalcs>);

/* Transmissive BCs, i.e. d/dx = 0 */
void TransBCs(Vec2D& u, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr);

/* Periodic BCs */
void PeriodicBCs(Vec2D& u, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr);

#endif
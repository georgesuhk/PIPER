#ifndef BC_HPP
#define BC_HPP

#include "sysCalcs.hpp"
#include "mesh.hpp"

// externally defined initial state, used in some BCs which hold one side fixed
extern Vec1D uLeftInit;


typedef void (*BCFunc)(Vec2D&, Mesh2D&, shared_ptr<SysCalcs>);

/* Transmissive BCs, i.e. d/dx = 0 */
void TransBCs(Vec2D& u, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr);

/* Periodic BCs */
void PeriodicBCs(Vec2D& u, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr);

/* Bohm BCs, used for detachment test */
void BohmBCs(Vec2D& u, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr);

void BohmBCs2(Vec2D& u, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr);

#endif
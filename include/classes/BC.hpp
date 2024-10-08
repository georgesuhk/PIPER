#ifndef BC_HPP
#define BC_HPP

#include "sysCalcs.hpp"
#include "mesh.hpp"

// externally defined initial state, used in some BCs which hold one side fixed
extern Vec1D uLeftInit;
extern Vec1D uRightInit;


typedef void (*BCFunc)(Vec2D&, Mesh2D&, shared_ptr<SysCalcs>);

/* Transmissive BCs, i.e. d/dx = 0 */
void TransBCs(Vec2D& u, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr);

/* Periodic BCs */
void PeriodicBCs(Vec2D& u, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr);

/* left reflective right transmissive*/
void LR_RT_BCs(Vec2D& u, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr);

/* top transmissive bottom reflective */
void TT_BR_BCs(Vec2D& u, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr);


/* X transmissive Y periodic BCs */
void XT_YP_BCs(Vec2D& u, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr);

/* Bohm BCs, used for detachment test */
void BohmBCs(Vec2D& u, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr);

void BohmBCs2(Vec2D& u, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr);

#endif
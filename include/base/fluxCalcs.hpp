#ifndef FLUXCALCS_HPP
#define FLUXCALCS_HPP

#include "mesh.hpp"
#include "sysCalcs.hpp"
#include "utils.hpp"

// following functions returns an array of numerical fluxes for the given scheme

/* Lax Friedrichs */
Vec2D getFluxLF(Vec2D& u, shared_ptr<SysCalcs> sysPtr, const Mesh2D& mesh, double dt, char axis);

/* Lax Friedrichs SLIC version */
Vec2D getFluxLF(Vec2D& uLHS, Vec2D& uRHS, shared_ptr<SysCalcs> sysPtr, const Mesh2D& mesh, double dt, char axis);

/* Richtmeyer */
Vec2D getFluxRI(Vec2D& u, shared_ptr<SysCalcs> sysPtr, const Mesh2D& mesh, double dt, char axis);

/* Richtmeyer SLIC version */
Vec2D getFluxRI(Vec2D& uLHS, Vec2D& uRHS, shared_ptr<SysCalcs> sysPtr, const Mesh2D& mesh, double dt, char axis);

#endif
#ifndef RECONSTRUCT_HPP
#define RECONSTRUCT_HPP

#include "utils.hpp"
#include "mesh.hpp"
#include "sysCalcs.hpp"
#include "BC.hpp"
#include "limiter.hpp"

/* returns 'deltas', the averaged difference between neighbouring variables */
Vec2D getDeltas(const Vec2D& u, Mesh2D& mesh, double w, char axis);

/* returns Vec2Ds of uL and uR which are the linearly reconstructed variables defined at cell boundary */
array<Vec2D,2> slopeRecon(const Vec2D& u, Mesh2D& mesh, double w, char axis, shared_ptr<SysCalcs> sysPtr, Limiter limFunc, BCFunc BCFunc);

/* half step update for reconstructed variables */
array<Vec2D,2> reconHS(Vec2D& uL, Vec2D& uR, shared_ptr<SysCalcs> sysPtr, Mesh2D& mesh, double dt, char axis, BCFunc BCFunc);



#endif
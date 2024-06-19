#ifndef INIT_HPP
#define INIT_HPP

#include "utils.hpp"
#include "sysCalcs.hpp"
#include "mesh.hpp"
#include "BC.hpp"

/* initialises a planar configuration with internal interfaces */
Vec2D initPlanar(vector<CellVec> initCellVecs, vector<double> interfacePositions, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr, BCFunc BCFunc, char axis);

#endif
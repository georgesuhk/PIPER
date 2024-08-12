#ifndef ORZAGINIT_HPP
#define ORZAGINIT_HPP

#include "utils.hpp"
#include "sysCalcs.hpp"
#include "mesh.hpp"
#include "BC.hpp"

Vec2D initOrzag(vector<double>& initParams, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr, BCFunc BCFunc);


#endif
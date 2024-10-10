#ifndef INIT2DTESTS_HPP
#define INIT2DTESTS_HPP

#include "utils.hpp"
#include "sysCalcs.hpp"
#include "mesh.hpp"
#include "BC.hpp"

Vec2D initOrzag(vector<double>& initParams, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr, BCFunc BCFunc);

Vec2D init_current_sheet(vector<double>& initParams, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr, BCFunc BCFunc);

Vec2D init_reconnection(vector<double>& initParams, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr, BCFunc BCFunc);

Vec2D init_mag_recon2(vector<double>& initParams, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr, BCFunc BCFunc);


#endif
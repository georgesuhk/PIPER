#ifndef DETACHINIT_HPP
#define DETACHINIT_HPP

#include "mutation++.h"
#include "utils.hpp"
#include "BC.hpp"
using namespace Mutation;
using namespace Mutation::Thermodynamics;
using namespace Mutation::Transport; 
using namespace std;

Vec2D initDetachment(vector<double>& initParams, string& mixName, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr, BCFunc BCFunc);



#endif

#ifndef SOURCEFUNCEX_HPP
#define SOURCEFUNCEX_HPP

#include "utils.hpp"
#include "sysCalcs.hpp"
#include "BC.hpp"
#include "calcs.hpp"

typedef Vec2D (*SourceFuncEx)(Vec2D&, Mesh2D&, shared_ptr<SysCalcs>, BCFunc);

/* relative ion - neutral velocity evolution */
Vec2D w_evolution_func(Vec2D& u, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr, BCFunc BC);


#endif
#include "utils.hpp"
#include "sysCalcs.hpp"

typedef void (*Exporter)(vector<Vec2D>&, vector<double>&, Mesh2D, shared_ptr<SysCalcs>, string);


/** 
 * Detailed exporter for fully ionised plasma
 * Contains:
 * Cell specific: rho, v, p, e, T, Bx, By, Bz
 * Total: U
*/
void ExportFIP(vector<Vec2D>& uTimeSeries, vector<double>& recordedTimes, Mesh2D mesh, shared_ptr<SysCalcs> sysPtr, string folder);

/** 
 * Detailed exporter for partially ionised plasma
 * Contains:
 * Cell specific: rho, v, p, e, T, Bx, By, Bz, wx, wy, wz, mass_frac_n, mass_frac_i
 * Total: U
*/
void ExportPIP(vector<Vec2D>& uTimeSeries, vector<double>& recordedTimes, Mesh2D mesh, shared_ptr<SysCalcs> sysPtr, string folder);
#include "BC.hpp"

void TransBCs(Vec2D& u, Mesh2D& mesh, shared_ptr<SysCalcs> SysPtr){
    
    //settings BCs for each inner vector (y direction)
    for (Vec1D& vec : u){
        vec[0] = vec[1];
        vec[mesh.nCellsY + 1] = vec[mesh.nCellsY];
    }

    //setting BCs for outer vector (x direction)
    u[0] = u[1];
    u[mesh.nCellsX + 1] = u[mesh.nCellsX];
}

void PeriodicBCs(Vec2D& u, Mesh2D& mesh, shared_ptr<SysCalcs> SysPtr){

    //settings BCs for each inner vector (y direction)
    for (Vec1D& vec : u){
        vec[0] = vec[mesh.nCellsY];
        vec[mesh.nCellsY + 1] = vec[1];
    }

    //setting BCs for outer vector (x direction)
    u[0] = u[mesh.nCellsX];
    u[mesh.nCellsX + 1] = u[1];
}
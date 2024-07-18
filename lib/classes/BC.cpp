#include "BC.hpp"

void TransBCs(Vec2D& u, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr){
    
    //settings BCs for each inner vector (y direction)
    for (Vec1D& vec : u){
        vec[0] = vec[1];
        vec[mesh.nCellsY + 1] = vec[mesh.nCellsY];
    }

    //setting BCs for outer vector (x direction)
    u[0] = u[1];
    u[mesh.nCellsX + 1] = u[mesh.nCellsX];
}

void LR_RT_BCs(Vec2D& u, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr){
    
    //settings BCs for each inner vector (y direction)
    for (Vec1D& vec : u){
        vec[0] = vec[1];
        vec[mesh.nCellsY + 1] = vec[mesh.nCellsY];
    }

    //setting BCs for outer vector (x direction)

    // left (reflective) ------
    u[0] = u[1];

    for (int j = 0; j < mesh.nCellsY+1; j++){
        u[0][j][1] = -u[1][j][1];
        u[0][j][6] = -u[1][j][6];
        u[0][j][9] = -u[1][j][9];
    }

    // right (transmissive)
    u[mesh.nCellsX + 1] = u[mesh.nCellsX];

}


void PeriodicBCs(Vec2D& u, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr){

    //settings BCs for each inner vector (y direction)
    for (Vec1D& vec : u){
        vec[0] = vec[mesh.nCellsY];
        vec[mesh.nCellsY + 1] = vec[1];
    }

    //setting BCs for outer vector (x direction)
    u[0] = u[mesh.nCellsX];
    u[mesh.nCellsX + 1] = u[1];
}


void BohmBCs(Vec2D& u, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr){

    // X DIRECTION ======
    double rho, vx, vy, vz, E0, e, p;

    /**
     * LEFT BOUNDARY
     * set to be fixed, based on value from initialisation
     * resembles a constant upstream heat flux
     */
    
    // for (int j = 0; j < mesh.nCellsY+2; j++){
    //     CellVec uPrim(cellVarsNums, 0);

    //     // information from last cell

    //     /* density */
    //     rho = u[1][j][0];
    //     vx = 0;
    //     vy = 0;
    //     vz = 0;

    //     /* pressure (power balance)*/
    //     E0 = sysPtr->get_E(uLeftInit[j]);
    //     e = E0 / rho;
    //     p = sysPtr->getEoSPtr()->interp_p(rho, e);

    //     // assembling

    //     uPrim[0] = rho;
    //     uPrim[1] = vx;
    //     uPrim[2] = vy;
    //     uPrim[3] = vz;
    //     uPrim[4] = p;
    //     uPrim[5] = u[1][j][5];
    //     uPrim[6] = u[1][j][6];
    //     uPrim[7] = u[1][j][7];
    //     uPrim[8] = u[1][j][8];
    //     // w BCs
    //     uPrim[9] = 0;
    //     uPrim[10] = 0;
    //     uPrim[11] = 0;
    //     u[0][j] = sysPtr->primToConserv(uPrim);
    // }
    
    // override to set left BC as constant
    u[0] = uLeftInit;



    /**
     * RIGHT BOUNDARY
     * sound speed outflow for vx, internal energy kept constant, and wx chosen to match flux
     */

    double Cs, mass_frac_n, mass_frac_i, n_n, n_i, wx;
    for (int j = 0; j < mesh.nCellsY+2; j++){
        CellVec uPrim(cellVarsNums, 0);

        // information from last cell

        /* density */
        rho = u[mesh.nCellsX][j][0];

        /* mass fractions */
        mass_frac_n = sysPtr->interp_mass_frac_n(u[mesh.nCellsX][j]);
        mass_frac_i = sysPtr->interp_mass_frac_i(u[mesh.nCellsX][j]);

        /* number densities */
        n_n = sysPtr->get_n_n(rho, mass_frac_n);
        n_i = sysPtr->get_n_i(rho, mass_frac_i);

        // creating BC values

        /* sound speed (for v_i via Bohm condition )*/
        Cs = sysPtr->get_Cs(u[mesh.nCellsX][j]);

        /* wx (particle balance) */
        wx = Cs * (1 + n_i / n_n);
        vx = Cs - mass_frac_n * wx;
        vy = 0;
        vz = 0;

        /* pressure (power balance) */
        E0 = sysPtr->get_E(uLeftInit[j]);
        e = (E0 - sysPtr->get_KE(rho, vx, vy, vz)) / rho;
        p = sysPtr->getEoSPtr()->interp_p(rho, e);

        // assembling

        uPrim[0] = rho;
        uPrim[1] = vx;
        uPrim[2] = 0;
        uPrim[3] = 0;
        uPrim[4] = p;
        
        // uPrim[4] = sysPtr->get_p(u[mesh.nCellsX][j]);
        uPrim[5] = u[mesh.nCellsX][j][5];
        uPrim[6] = u[mesh.nCellsX][j][6];
        uPrim[7] = u[mesh.nCellsX][j][7];
        uPrim[8] = u[mesh.nCellsX][j][8];
        // w BCs
        uPrim[9] = wx;
        uPrim[10] = 0;
        uPrim[11] = 0;
        u[mesh.nCellsX+1][j] = sysPtr->primToConserv(uPrim);
    }



    // Y DIRECTION ======
    for (Vec1D& vec : u){
        vec[0] = vec[1];
        vec[mesh.nCellsY + 1] = vec[mesh.nCellsY];
    }
}



void BohmBCs2(Vec2D& u, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr){

    // X DIRECTION ======
    double rho, vx, vy, vz, E0, e, p;

    /**
     * LEFT BOUNDARY
     * set to be fixed, based on value from initialisation
     * resembles a constant upstream heat flux
     */
    
    for (int j = 0; j < mesh.nCellsY+2; j++){
        CellVec uPrim(cellVarsNums, 0);

        // information from last cell

        /* density */
        rho = u[1][j][0];
        vx = u[1][j][1] / rho;
        vy = u[1][j][2] / rho;
        vz = u[1][j][3] / rho;

        /* pressure (power balance)*/
        p = sysPtr->get_p(u[1][j]);

        // assembling

        uPrim[0] = rho;
        uPrim[1] = -vx;
        uPrim[2] = vy;
        uPrim[3] = vz;
        uPrim[4] = p;
        uPrim[5] = u[1][j][5];
        uPrim[6] = u[1][j][6];
        uPrim[7] = u[1][j][7];
        uPrim[8] = u[1][j][8];
        // w BCs
        uPrim[9] = u[1][j][9];
        uPrim[10] = u[1][j][10];
        uPrim[11] = u[1][j][11];
        u[0][j] = sysPtr->primToConserv(uPrim);
    }
    
    // override to set left BC as constant
    // u[0] = uLeftInit;



    /**
     * RIGHT BOUNDARY
     * sound speed outflow for vx, internal energy kept constant, and wx chosen to match flux
     */

    double Cs, mass_frac_n, mass_frac_i, n_n, n_i, wx;
    for (int j = 0; j < mesh.nCellsY+2; j++){
        CellVec uPrim(cellVarsNums, 0);

        // information from last cell

        /* density */
        rho = u[mesh.nCellsX][j][0];

        /* mass fractions */
        mass_frac_n = sysPtr->interp_mass_frac_n(u[mesh.nCellsX][j]);
        mass_frac_i = sysPtr->interp_mass_frac_i(u[mesh.nCellsX][j]);

        /* number densities */
        n_n = sysPtr->get_n_n(rho, mass_frac_n);
        n_i = sysPtr->get_n_i(rho, mass_frac_i);

        // creating BC values

        /* sound speed (for v_i via Bohm condition )*/
        Cs = sysPtr->get_Cs(u[mesh.nCellsX][j]);

        /* wx (particle balance) */
        wx = Cs * (1 + n_i / n_n);
        vx = Cs - mass_frac_n * wx;
        vy = 0;
        vz = 0;

        /* pressure (power balance) */
        E0 = sysPtr->get_E(uLeftInit[j]);
        e = (E0 - sysPtr->get_KE(rho, vx, vy, vz)) / rho;
        p = sysPtr->getEoSPtr()->interp_p(rho, e);

        // assembling

        uPrim[0] = rho;
        uPrim[1] = vx;
        uPrim[2] = 0;
        uPrim[3] = 0;
        // uPrim[4] = p;
        uPrim[4] = sysPtr->get_p(u[mesh.nCellsX][j]);
        uPrim[5] = u[mesh.nCellsX][j][5];
        uPrim[6] = u[mesh.nCellsX][j][6];
        uPrim[7] = u[mesh.nCellsX][j][7];
        uPrim[8] = u[mesh.nCellsX][j][8];
        // w BCs
        uPrim[9] = wx;
        uPrim[10] = 0;
        uPrim[11] = 0;
        u[mesh.nCellsX+1][j] = sysPtr->primToConserv(uPrim);
    }



    // Y DIRECTION ======
    for (Vec1D& vec : u){
        vec[0] = vec[1];
        vec[mesh.nCellsY + 1] = vec[mesh.nCellsY];
    }
}
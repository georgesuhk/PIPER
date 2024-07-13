#include "sourceFuncEx.hpp"


Vec2D heating(Vec2D& u, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr, BCFunc BC){
    Vec2D S = makeVec2D(mesh.nCellsX + 2, mesh.nCellsY + 2);

    for (int i = 1; i < round(mesh.nCellsX/2); i++){
        for (int j = 1; j < mesh.nCellsY+1; j++){
            CellVec cellVec(12, 0);
            cellVec[4] = 80000;

            S[i][j] = cellVec;
        }
    }

    return S;
}


Vec2D w_evolution_func(Vec2D& u, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr, BCFunc BC){
    Vec2D S = makeVec2D(mesh.nCellsX + 2, mesh.nCellsY + 2);

    double rho, vx, vy, vz, wx, wy, wz;
    double dwx_dx, dwx_dy, dwy_dx, dwy_dy, dwz_dx, dwz_dy;
    double dvx_dx, dvx_dy, dvy_dx, dvy_dy, dvz_dx, dvz_dy;
    double v_dot_nabla_w_x, v_dot_nabla_w_y, v_dot_nabla_w_z;
    double w_dot_nabla_v_x, w_dot_nabla_v_y, w_dot_nabla_v_z;
    double w_dot_nabla_w_x, w_dot_nabla_w_y, w_dot_nabla_w_z;
    double mass_frac_n, mass_frac_i, n_n, n_i, m_n, m_i, alpha_in, alpha_en, T, prefactor;
    for (int i = 1; i < mesh.nCellsX+1; i++){
        for (int j = 1; j < mesh.nCellsY+1; j++){
            CellVec cellVec(12, 0);

            rho = u[i][j][0];
            vx = u[i][j][1] / rho;
            vy = u[i][j][2] / rho;
            vz = u[i][j][3] / rho;
            wx = u[i][j][9];
            wy = u[i][j][10];
            wz = u[i][j][11];

            /* components of the convective derivative */

            /* wx derivatives */
            dwx_dx = ( u[i+1][j][9] - u[i-1][j][9] ) / (2*mesh.dx);
            dwx_dy = ( u[i][j+1][9] - u[i][j-1][9] ) / (2*mesh.dy);

            /* wy derivatives */
            dwy_dx = ( u[i+1][j][10] - u[i-1][j][10] ) / (2*mesh.dx);
            dwy_dy = ( u[i][j+1][10] - u[i][j-1][10] ) / (2*mesh.dy);

            /* wz derivatives */
            dwz_dx = ( u[i+1][j][11] - u[i-1][j][11] ) / (2*mesh.dx);
            dwz_dy = ( u[i][j+1][11] - u[i][j-1][11] ) / (2*mesh.dy);

            /* vx derivatives */
            dvx_dx = ( u[i+1][j][1] / u[i+1][j][0] - u[i-1][j][1] / u[i-1][j][0])  / (2*mesh.dx);
            dvx_dy = ( u[i][j+1][1] / u[i][j+1][0] - u[i][j-1][1] / u[i][j-1][0])  / (2*mesh.dy);

            /* vy derivatives */
            dvy_dx = ( u[i+1][j][2] / u[i+1][j][0] - u[i-1][j][2] / u[i-1][j][0])  / (2*mesh.dx);
            dvy_dy = ( u[i][j+1][2] / u[i][j+1][0] - u[i][j-1][2] / u[i][j-1][0])  / (2*mesh.dy);

            /* vz derivatives */
            dvz_dx = ( u[i+1][j][3] / u[i+1][j][0] - u[i-1][j][3] / u[i-1][j][0])  / (2*mesh.dx);
            dvx_dy = ( u[i][j+1][3] / u[i][j+1][0] - u[i][j-1][3] / u[i][j-1][0])  / (2*mesh.dy);

            /* compound terms */
            v_dot_nabla_w_x = vx * dwx_dx + vy * dwx_dy;
            v_dot_nabla_w_y = vx * dwy_dx + vy * dwy_dy;
            v_dot_nabla_w_z = vx * dwz_dx + vy * dwz_dy;

            w_dot_nabla_v_x = wx * dvx_dx + wy * dvx_dy;
            w_dot_nabla_v_y = wx * dvy_dx + wy * dvy_dy;
            w_dot_nabla_v_z = wx * dvz_dx + wy * dvz_dy;

            w_dot_nabla_w_x = wx * dwx_dx + wy * dwx_dy;
            w_dot_nabla_w_y = wx * dwy_dx + wy * dwy_dy;
            w_dot_nabla_w_z = wx * dwz_dx + wy * dwz_dy;

            // collision coefficients 
            // mass_frac_i = sysPtr->interp_mass_frac_i(u[i][j]);
            mass_frac_n = sysPtr->interp_mass_frac_n(u[i][j]);
            // n_i = sysPtr->get_n_i(rho, mass_frac_i);
            // n_n = sysPtr->get_n_n(rho, mass_frac_n); 
            // m_n = sysPtr->get_m_n();
            // m_i = sysPtr->get_m_i();
            // T = sysPtr->get_T(u[i][j]);
            // alpha_in = get_coll_coeff_in(n_i, n_n, m_i, m_n, T);
            // alpha_en = get_coll_coeff_en(n_i, n_n, m_i, m_n, T);

            /* size comparison */
            if (i == round(mesh.nCellsX/2) && j == 1){
                cout << endl;
                cout << "v_dot_nabla_w: " << v_dot_nabla_w_x << endl;
                cout << "w_dot_nabla_v: " << w_dot_nabla_v_x << endl;
                cout << "w_dot_nabla_w: " << w_dot_nabla_w_x << endl;
                cout << "w_dot_nabla_w * mass_frac_n: " << mass_frac_n * w_dot_nabla_w_x << endl;
            }


            // assembling

            // temp test without second evolutionary term
            cellVec[9] = - v_dot_nabla_w_x - w_dot_nabla_v_x - mass_frac_n * w_dot_nabla_w_x;
            cellVec[10] = - v_dot_nabla_w_y - w_dot_nabla_v_y - mass_frac_n * w_dot_nabla_w_y;
            cellVec[11] = - v_dot_nabla_w_z - w_dot_nabla_v_z - mass_frac_n * w_dot_nabla_w_z;

            S[i][j] = cellVec;
        }
    }

    // cout << "S: " << endl;
    // cout << S[round(mesh.nCellsX / 2)][1] << endl;

    return S;
}

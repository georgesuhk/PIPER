#include "sourceFuncEx.hpp"


Vec2D heating(Vec2D& u, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr, BCFunc BC){
    Vec2D S = makeVec2D(mesh.nCellsX + 2, mesh.nCellsY + 2);

    for (int i = 1; i < round(mesh.nCellsX/2); i++){
        for (int j = 1; j < mesh.nCellsY+1; j++){
            CellVec cellVec(12, 0);
            cellVec[4] = 50000;

            S[i][j] = cellVec;
        }
    }

    return S;
}


Vec2D w_evolution_func(Vec2D& u, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr, BCFunc BC){
    Vec2D S = makeVec2D(mesh.nCellsX + 2, mesh.nCellsY + 2);

    /* constant values */
    double two_dx = 2 * mesh.dx;
    double two_dy = 2 * mesh.dy;

    /* state vars */
    double rho, vx, vy, vz, Bx, By, Bz, wx, wy, wz;
    double T, n_n, n_i, m_n, m_i;
    double mass_frac_i, mass_frac_n;

    /* coll integrals */
    double alpha_en, alpha_in, coll_term_prefac;

    /* state vars used in derivatives */
    double mfi_x_plus, mfi_x_minus, mfi_y_plus, mfi_y_minus;
    double T_x_plus, T_x_minus, T_y_plus, T_y_minus;

    /* derivatives */
    double dwx_dx, dwx_dy, dwy_dx, dwy_dy, dwz_dx, dwz_dy;
    double dvx_dx, dvx_dy, dvy_dx, dvy_dy, dvz_dx, dvz_dy;
    double dBx_dx, dBx_dy, dBy_dx, dBy_dy, dBz_dx, dBz_dy;
    double d_mfi_dx, d_mfi_dy;
    double d_rho_i_dx, d_rho_i_dy;
    double d_pn_dx, d_pn_dy, d_pei_dx, d_pei_dy;

    /* compound terms */
    double v_dot_nabla_w_x, v_dot_nabla_w_y, v_dot_nabla_w_z;
    double w_dot_nabla_v_x, w_dot_nabla_v_y, w_dot_nabla_v_z;
    double w_dot_nabla_w_x, w_dot_nabla_w_y, w_dot_nabla_w_z;
    double w_dot_d_mfi;
    double mag_prefac, mag_term_x, mag_term_y, mag_term_z;
    double p_prefac, p_term_x, p_term_y, p_term_z;


    for (int i = 1; i < mesh.nCellsX+1; i++){
        for (int j = 1; j < mesh.nCellsY+1; j++){
            CellVec cellVec(12, 0);

            rho = u[i][j][0];
            vx = u[i][j][1] / rho;
            vy = u[i][j][2] / rho;
            vz = u[i][j][3] / rho;
            Bx = u[i][j][5];
            By = u[i][j][6];
            Bz = u[i][j][7];
            wx = u[i][j][9];
            wy = u[i][j][10];
            wz = u[i][j][11];

            T = sysPtr->get_T(i, j);

            T_x_plus = sysPtr->get_T(i+1, j);
            T_x_minus = sysPtr->get_T(i-1, j);
            T_y_plus = sysPtr->get_T(i, j+1);
            T_y_minus = sysPtr->get_T(i, j-1);

            mass_frac_i = sysPtr->get_mass_frac_i(i, j);
            mass_frac_n = 1 - mass_frac_i;
            mfi_x_plus = sysPtr->get_mass_frac_i(i+1, j);
            mfi_x_minus = sysPtr->get_mass_frac_i(i-1, j);
            mfi_y_plus = sysPtr->get_mass_frac_i(i, j+1);
            mfi_y_minus = sysPtr->get_mass_frac_i(i, j-1);


            /* components of the convective derivative */

            /* wx derivatives */
            dwx_dx = ( u[i+1][j][9] - u[i-1][j][9] ) / two_dx;
            dwx_dy = ( u[i][j+1][9] - u[i][j-1][9] ) / two_dy;

            /* wy derivatives */
            dwy_dx = ( u[i+1][j][10] - u[i-1][j][10] ) / two_dx;
            dwy_dy = ( u[i][j+1][10] - u[i][j-1][10] ) / two_dy;

            /* wz derivatives */
            dwz_dx = ( u[i+1][j][11] - u[i-1][j][11] ) / two_dx;
            dwz_dy = ( u[i][j+1][11] - u[i][j-1][11] ) / two_dy;

            /* vx derivatives */
            dvx_dx = ( u[i+1][j][1] / u[i+1][j][0] - u[i-1][j][1] / u[i-1][j][0])  / two_dx;
            dvx_dy = ( u[i][j+1][1] / u[i][j+1][0] - u[i][j-1][1] / u[i][j-1][0])  / two_dy;

            /* vy derivatives */
            dvy_dx = ( u[i+1][j][2] / u[i+1][j][0] - u[i-1][j][2] / u[i-1][j][0])  / two_dx;
            dvy_dy = ( u[i][j+1][2] / u[i][j+1][0] - u[i][j-1][2] / u[i][j-1][0])  / two_dy;

            /* vz derivatives */
            dvz_dx = ( u[i+1][j][3] / u[i+1][j][0] - u[i-1][j][3] / u[i-1][j][0])  / two_dx;
            dvx_dy = ( u[i][j+1][3] / u[i][j+1][0] - u[i][j-1][3] / u[i][j-1][0])  / two_dy;

            /* Bx derivatives */
            dBx_dx = ( u[i+1][j][5] - u[i-1][j][5] ) / two_dx;
            dBx_dy = ( u[i][j+1][5] - u[i][j-1][5] ) / two_dy;

            /* Bx derivatives */
            dBy_dx = ( u[i+1][j][6] - u[i-1][j][6] ) / two_dx;
            dBy_dy = ( u[i][j+1][6] - u[i][j-1][6] ) / two_dy;

            /* Bz derivatives */
            dBz_dx = ( u[i+1][j][7] - u[i-1][j][7] ) / two_dx;
            dBz_dy = ( u[i][j+1][7] - u[i][j-1][7] ) / two_dy;

            /* mass_frac_i derivatives */
            d_mfi_dx = ( mfi_x_plus - mfi_x_minus ) / two_dx;
            d_mfi_dy = ( mfi_y_plus - mfi_y_minus ) / two_dy;

            d_rho_i_dx = ( mfi_x_plus * u[i+1][j][0] - mfi_x_minus * u[i-1][j][0] ) / two_dx;
            d_rho_i_dy = ( mfi_y_plus * u[i][j+1][0] - mfi_y_minus * u[i][j-1][0] ) / two_dy;

            /* partial pressure derivatives */
            d_pn_dx = kBScaled * ( sysPtr->get_n_n(u[i+1][j][0], 1-mfi_x_plus) * T_x_plus - sysPtr->get_n_n(u[i-1][j][0], 1-mfi_x_minus) * T_x_minus) / two_dx;
            d_pn_dy = kBScaled * ( sysPtr->get_n_n(u[i][j+1][0], 1-mfi_y_plus) * T_y_plus - sysPtr->get_n_n(u[i][j-1][0], 1-mfi_y_minus) * T_y_minus) / two_dy;

            d_pei_dx = kBScaled * ( 2 * sysPtr->get_n_i(u[i+1][j][0], mfi_x_plus) * T_x_plus - 2 * sysPtr->get_n_i(u[i-1][j][0], mfi_x_minus) * T_x_minus) / two_dx;
            d_pei_dy = kBScaled * ( 2 * sysPtr->get_n_i(u[i][j+1][0], mfi_y_plus) * T_y_plus - 2 * sysPtr->get_n_i(u[i][j-1][0], mfi_y_minus) * T_y_minus) / two_dy;
            

            /* compound terms */

            // term 1
            v_dot_nabla_w_x = vx * dwx_dx + vy * dwx_dy;
            v_dot_nabla_w_y = vx * dwy_dx + vy * dwy_dy;
            v_dot_nabla_w_z = vx * dwz_dx + vy * dwz_dy;

            // term 2
            w_dot_nabla_v_x = wx * dvx_dx + wy * dvx_dy;
            w_dot_nabla_v_y = wx * dvy_dx + wy * dvy_dy;
            w_dot_nabla_v_z = wx * dvz_dx + wy * dvz_dy;

            // term 3
            w_dot_nabla_w_x = wx * dwx_dx + wy * dwx_dy;
            w_dot_nabla_w_y = wx * dwy_dx + wy * dwy_dy;
            w_dot_nabla_w_z = wx * dwz_dx + wy * dwz_dy;

            // term 4
            w_dot_d_mfi = wx * d_mfi_dx + wy * d_mfi_dy;

            // term 5 (magnetic term)
            mag_prefac = 1 / (mass_frac_i * rho + 1e-200);
            mag_term_x = mag_prefac * (By * dBx_dy - By * dBy_dx - Bz * dBz_dx);
            mag_term_y = mag_prefac * (Bx * dBy_dx - Bx * dBx_dy - Bz * dBz_dy);
            mag_term_z = mag_prefac * (Bx * dBz_dx + By * dBz_dy);

            // term 6 (p term)
            p_prefac = 1 / (mass_frac_i * mass_frac_n * rho + 1e-200);
            p_term_x = p_prefac * (mass_frac_i * d_pn_dx - mass_frac_n * d_pei_dx);
            p_term_y = p_prefac * (mass_frac_i * d_pn_dy - mass_frac_n * d_pei_dy);
            p_term_z = 0;


            // collision coefficients 
            n_i = sysPtr->get_n_i(rho, mass_frac_i);
            n_n = sysPtr->get_n_n(rho, mass_frac_n); 
            m_n = sysPtr->get_m_n();
            m_i = sysPtr->get_m_i();
            alpha_in = get_coll_coeff_in(n_i, n_n, m_i, m_n, T);
            alpha_en = get_coll_coeff_en(n_i, n_n, m_i, m_n, T);
            coll_term_prefac = (alpha_in + alpha_en) / (rho * mass_frac_i * mass_frac_n);


            /* size comparison */
            // if (i == 10 && j == 1){
            //     cout << endl;
            //     cout << "v_dot_nabla_w: " << -v_dot_nabla_w_x << endl;
            //     cout << "w_dot_nabla_v: " << -w_dot_nabla_v_x << endl;
            //     cout << "w_dot_nabla_w: " << -w_dot_nabla_w_x << endl;
            //     cout << "(1 - 2*mass_frac_i) * w dot nabla w: " << -(1 - 2*mass_frac_i) * w_dot_nabla_w_x << endl;
            //     cout << "dmfi term: " << wx * w_dot_d_mfi << endl;
            //     cout << "mag term: " << mag_term_x << endl;
            //     cout << "p term: " << p_term_x << endl;
            //     cout << "coll term: " << -coll_term_prefac * wx << endl; 
            // }


            // assembling
            double dampFac = coll_term_prefac * 0.0;

            // cellVec[9] = - v_dot_nabla_w_x - w_dot_nabla_v_x - (1 - 2*mass_frac_i) * w_dot_nabla_w_x + wx * w_dot_d_mfi + mag_term_x + p_term_x - dampFac * wx;
            // cellVec[10] = - v_dot_nabla_w_y - w_dot_nabla_v_y - (1 - 2*mass_frac_i) * w_dot_nabla_w_y + wy * w_dot_d_mfi + mag_term_y + p_term_y - dampFac * wy;
            // cellVec[11] = - v_dot_nabla_w_z - w_dot_nabla_v_z - (1 - 2*mass_frac_i) * w_dot_nabla_w_z + wz * w_dot_d_mfi + mag_term_z + p_term_z - dampFac * wz;

            cellVec[9] = -v_dot_nabla_w_x - w_dot_nabla_v_x - (1 - 2*mass_frac_i) * w_dot_nabla_w_x + wx * w_dot_d_mfi + mag_term_x + 1*p_term_x - dampFac * wx;
            cellVec[10] = -v_dot_nabla_w_y - w_dot_nabla_v_y - (1 - 2*mass_frac_i) * w_dot_nabla_w_y + wy * w_dot_d_mfi + mag_term_y + 1*p_term_y - dampFac * wy;
            cellVec[11] = -v_dot_nabla_w_z - w_dot_nabla_v_z - (1 - 2*mass_frac_i) * w_dot_nabla_w_z + wz * w_dot_d_mfi + mag_term_z + 1*p_term_z - dampFac * wz;

            // if (i == round(mesh.nCellsX/2) && j == 1){
            //     cout << cellVec[9] << endl;
            // }

            S[i][j] = cellVec;
        }
    }

    // cout << "S: " << endl;
    // cout << S[round(mesh.nCellsX / 2)][1] << endl;

    return S;
}

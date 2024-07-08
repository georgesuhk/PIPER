#include "sourceFuncEx.hpp"


Vec2D w_evolution_func(Vec2D& u, Mesh2D& mesh, shared_ptr<SysCalcs> sysPtr, BCFunc BC){
    Vec2D S = makeVec2D(mesh.nCellsX + 2, mesh.nCellsY + 2);

    CellVec cellVec(12, 0);
    double rho, wx, wy, wz;
    double dwx_dx, dwx_dy, dwy_dx, dwy_dy, dwz_dx, dwz_dy;
    double mass_frac_n, mass_frac_i, n_n, n_i, m_n, m_i, alpha_in, alpha_en, T, prefactor;
    for (int i = 1; i < mesh.nCellsX+1; i++){
        for (int j = 1; j < mesh.nCellsY+1; j++){

            rho = u[i][j][0];

            //dw/dt source terms

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

            // collision coefficients 
            mass_frac_i = sysPtr->interp_mass_frac_i(u[i][j]);
            mass_frac_n = sysPtr->interp_mass_frac_n(u[i][j]);
            n_i = sysPtr->get_n_i(rho, mass_frac_i);
            n_n = sysPtr->get_n_n(rho, mass_frac_n); 
            m_n = sysPtr->get_m_n();
            m_i = sysPtr->get_m_i();
            T = sysPtr->get_T(u[i][j]);
            alpha_in = get_coll_coeff_in(n_i, n_n, m_i, m_n, T);
            alpha_en = get_coll_coeff_en(n_i, n_n, m_i, m_n, T);

            /* prefactor for second source term */
            prefactor = (alpha_in * alpha_en) / (rho * mass_frac_i * mass_frac_n);

            // if (i == 50 && j == 1){
            //     cout << "dwx_dx: " << dwx_dx << endl;
            //     cout << "dwx_dy: " << dwx_dy << endl;
            //     cout << u[i+1][j][9] << endl;
            //     cout << u[i+1][j][9] << endl;

            // }

            // cout << "dwy_dx: " << dwy_dx << endl;
            // cout << "dwy_dx: " << dwy_dy << endl;
            // cout << "dwz_dx: " << dwz_dx << endl;
            // cout << "dwz_dy: " << dwz_dy << endl;

            // assembling
            // cellVec[9] = wx * dwx_dx + wy * dwx_dy + prefactor * wx;
            // cellVec[10] = wy * dwy_dx + wy * dwy_dy + prefactor * wy;
            // cellVec[11] = wz * dwz_dx + wz * dwz_dy + prefactor * wz;

            // temp test without second evolutionary term
            cellVec[9] = wx * dwx_dx + wy * dwx_dy;
            cellVec[10] = wy * dwy_dx + wy * dwy_dy;
            cellVec[11] = wz * dwz_dx + wz * dwz_dy;

            S[i][j] = cellVec;
        }
    }

    cout << "S: " << endl;
    cout << S[50][1] << endl;

    BC(S, mesh, sysPtr);
    return S;
}

#include "implicitSolvers.hpp"

// DIFFUSION COEFFICIENTS ======
double get_Resis(CellVec& u, int i, int j, shared_ptr<SysCalcs> sysPtr, bool interp){
    return sysPtr->get_Resis(u, i, j, interp);
}



// IMPLICIT SOLVERS ======

/**
 * Crank Nicolson solver for the diffusion equation
*/

void CN_Diffusion_Solver(Vec2D& u, int varIdx, GetCoeffFunc getCoeffFunc, shared_ptr<SysCalcs> sysPtr, Mesh2D& mesh, double& dt, string diffBCType, BCFunc BCFunc){
    
    /* intermediate matix for after x update */
    Vec2D uMid = u;

    /* diffusion coefficient */
    double diffCoeff, diffCoeff_next, diffCoeff_prev;

    /* coefficent not containing diffCoeff*/
    double alpha_x = (dt/2.0) / (2 * pow(mesh.dx, 2.0)); 
    double alpha_y = (dt/2.0) / (2 * pow(mesh.dy, 2.0)); 

    /* whether diffusion coefficient should be reinterpreted */
    bool interp = false; 

    double b_vec_term;
    vector<double> solution;

    // X update ------

    // Loop over rows
    for (int j = 1; j < mesh.nCellsY+1; j++){
        
        /* matrix components */
        vector<double> A_diag_X, A_off_diag_l_X, A_off_diag_r_X, b_vec_X;

        // matrix construction
        for (int i = 1; i < mesh.nCellsX+1; i++){

            // obtaining cell specific coefficients
            diffCoeff = getCoeffFunc(u[i][j], i, j, sysPtr, interp);
            diffCoeff_next = getCoeffFunc(u[i][j+1], i, j+1, sysPtr, interp);
            diffCoeff_prev = getCoeffFunc(u[i][j-1], i, j-1, sysPtr, interp);

            // creating A matrix 
            A_diag_X.push_back(1 + 2 * diffCoeff * alpha_x);
            A_off_diag_l_X.push_back(-1 * diffCoeff * alpha_x);
            A_off_diag_r_X.push_back(-1 * diffCoeff * alpha_x);

            // creating B matrix
            if (j == 1){
                if (diffBCType == "Neumann"){
                    b_vec_term = (1 - (2/3) * alpha_y * diffCoeff) * u[i][j][varIdx] +  ((2/3) * alpha_y * diffCoeff_next * u[i][j+1][varIdx]); 
                } else {
                    throw runtime_error("diffBCType: " + diffBCType + " is not available.");
                }
            } else if (j == mesh.nCellsY){
                if (diffBCType == "Neumann"){
                    b_vec_term = (1 - (2/3) * alpha_y * diffCoeff) * u[i][j][varIdx] +  ((2/3) * alpha_y * diffCoeff_prev * u[i][j-1][varIdx]); 
                } else {
                    throw runtime_error("diffBCType: " + diffBCType + " is not available.");
                }                
            } else {
                b_vec_term = (1 - 2 * alpha_y * diffCoeff) * u[i][j][varIdx] + alpha_y * (diffCoeff_prev * u[i][j-1][varIdx] + diffCoeff_next * u[i][j+1][varIdx]);
            }

            b_vec_X.push_back(b_vec_term);
        }

        A_off_diag_l_X.pop_back();
        A_off_diag_r_X.pop_back(); 

        // Modifiying A matrix elements according to BCs

        if (diffBCType == "Neumann"){
            A_diag_X[0] = 1 + (2.0/3.0) * alpha_x * getCoeffFunc(u[1][j], 1, j, sysPtr, interp);
            A_diag_X.back() = 1 + (2.0/3.0)  * alpha_x * getCoeffFunc(u[mesh.nCellsX][j], mesh.nCellsX, j, sysPtr, interp);
            A_off_diag_r_X[0] = -1 * (2.0/3.0) * alpha_x * getCoeffFunc(u[2][j], 2, j, sysPtr, interp);
            A_off_diag_l_X.back() = -1 * (2.0/3.0) * alpha_x * getCoeffFunc(u[mesh.nCellsX-1][j], mesh.nCellsX-1, j, sysPtr, interp);
        } else {
            throw runtime_error("diffBCType: " + diffBCType + " is not available.");
        }

        /* solution vector, with u1 term at the zeroth index, hence need i-1 shift to match u */
        solution = thomasAlgo(mesh, A_diag_X, A_off_diag_l_X, A_off_diag_r_X, b_vec_X, 'x');

        // Applying solution to u
        for (int i = 1; i < mesh.nCellsX+1; i++){
            uMid[i][j][varIdx] = solution[i-1];
        }
    }

    // Y update ------
    interp = true;

    // Loop over columns
    for (int i = 1; i < mesh.nCellsX+1; i++){

        /* matrix components */
        vector<double> A_diag_Y, A_off_diag_l_Y, A_off_diag_r_Y, b_vec_Y;

        // Matrix construction
        for (int j = 1; j < mesh.nCellsY+1; j++){

            //obtaining cell specific coefficients
            diffCoeff = getCoeffFunc(uMid[i][j], i, j, sysPtr, interp);
            diffCoeff_next = getCoeffFunc(uMid[i+1][j], i+1, j, sysPtr, interp);
            diffCoeff_prev = getCoeffFunc(uMid[i-1][j], i-1, j, sysPtr, interp);

            // creating A matrix 
            A_diag_Y.push_back(1 + 2 * diffCoeff * alpha_y);
            A_off_diag_l_Y.push_back(-1 * diffCoeff * alpha_y);
            A_off_diag_r_Y.push_back(-1 * diffCoeff * alpha_y);

            // creating B matrix
            if (i == 1){
                if (diffBCType == "Neumann"){
                    b_vec_term = (1 - (2/3) * alpha_x * diffCoeff) * uMid[i][j][varIdx] +  ((2/3) * alpha_x * diffCoeff_next * uMid[i+1][j][varIdx]); 
                } else {
                    throw runtime_error("diffBCType: " + diffBCType + " is not available.");
                }
            } else if (i == mesh.nCellsX){
                if (diffBCType == "Neumann"){
                    b_vec_term = (1 - (2/3) * alpha_x * diffCoeff) * uMid[i][j][varIdx] +  ((2/3) * alpha_x * diffCoeff_prev * uMid[i-1][j][varIdx]); 
                } else {
                    throw runtime_error("diffBCType: " + diffBCType + " is not available.");
                }                
            } else {
                b_vec_term = (1 - 2 * alpha_x * diffCoeff) * uMid[i][j][varIdx] + alpha_x * (diffCoeff_prev * uMid[i-1][j][varIdx] + diffCoeff_next * uMid[i+1][j][varIdx]);
            }

            b_vec_Y.push_back(b_vec_term);
        }

        A_off_diag_l_Y.pop_back();
        A_off_diag_r_Y.pop_back(); 

        // Modifiying A matrix elements according to BCs
        if (diffBCType == "Neumann"){
            A_diag_Y[0] = 1 + (2.0/3.0) * alpha_y * getCoeffFunc(uMid[i][1], i, 1, sysPtr, interp);
            A_diag_Y.back() = 1 + (2.0/3.0)  * alpha_y * getCoeffFunc(uMid[i][mesh.nCellsY], i, mesh.nCellsY, sysPtr, interp);
            A_off_diag_r_Y[0] = -1 * (2.0/3.0) * alpha_y * getCoeffFunc(uMid[i][2], i, 2, sysPtr, interp);
            A_off_diag_l_Y.back() = -1 * (2.0/3.0) * alpha_y * getCoeffFunc(uMid[i][mesh.nCellsY-1], i, mesh.nCellsY-1, sysPtr, interp);
        } else {
            throw runtime_error("diffBCType: " + diffBCType + " is not available.");
        }

        /* solution vector, with u1 term at the zeroth index, hence need i-1 shift to match u */
        solution = thomasAlgo(mesh, A_diag_Y, A_off_diag_l_Y, A_off_diag_r_Y, b_vec_Y, 'y');

        // Applying solution to u
        for (int j = 1; j < mesh.nCellsY+1; j++){
            u[i][j][varIdx] = solution[j-1];
        }
    }

    BCFunc(u, mesh, sysPtr);
    sysPtr->getEoSPtr()->cacheAll(u, mesh);
}




#include "implicitSolvers.hpp"

// DIFFUSION COEFFICIENTS ======
double get_Resis(CellVec& u, int i, int j, shared_ptr<SysCalcs> sysPtr, bool interp){
    return sysPtr->get_Resis(u, i, j, interp);
}

double get_therm_con(CellVec& u, int i, int j, shared_ptr<SysCalcs> sysPtr, bool interp){
    return sysPtr->interp_therm_con(u);
}



// IMPLICIT SOLVERS ======

/**
 * Crank Nicolson solver for the diffusion equation, applied to the resistive source term
*/

void CN_Diffusion_Solver_ADI(Vec2D& u, vector<int> varIdxList, GetCoeffFunc getCoeffFunc, shared_ptr<SysCalcs> sysPtr, Mesh2D& mesh, double& dt, string diffBCType, BCFunc BCFunc){
    
    /* intermediate matix for after x update */
    Vec2D uMid = u;
    double varIdx = varIdxList[0];

    /* diffusion coefficient */
    double diffCoeff, diffCoeff_next, diffCoeff_prev;

    /* coefficent not containing diffCoeff*/
    double alpha_x = (dt/2.0) / (2 * pow(mesh.dx, 2.0)); 
    double alpha_y = (dt/2.0) / (2 * pow(mesh.dy, 2.0)); 

    /* whether diffusion coefficient should be reinterpreted */
    bool interp = true; 

    double b_vec_term;
    vector<double> solution;

    // Y update ------

    // Loop over rows
    for (int j = 1; j < mesh.nCellsY+1; j++){
        
        /* matrix components */
        vector<double> A_diag_X, A_off_diag_l_X, A_off_diag_r_X, b_vec_X;

        // matrix construction
        for (int i = 1; i < mesh.nCellsX+1; i++){

            // obtaining cell specific coefficients
            // diffCoeff = getCoeffFunc(u[i][j], i, j, sysPtr, interp);
            // diffCoeff_next = getCoeffFunc(u[i][j+1], i, j+1, sysPtr, interp);
            // diffCoeff_prev = getCoeffFunc(u[i][j-1], i, j-1, sysPtr, interp);
            double Lr = 0.05;
            double x = indexToPos(mesh, i, 'x');
            double y = indexToPos(mesh, j, 'y');
            double diffCoeffInterest = 1.0 * (cos(myPI * x / (2*Lr))+1) * (cos(myPI*y/(4 * Lr))+1)/4;;
            if (fabs(x) < Lr && fabs(y) < 2 * Lr){
                diffCoeff = diffCoeffInterest;
                diffCoeff_next = diffCoeffInterest;
                diffCoeff_prev = diffCoeffInterest;
            } else {
                diffCoeff = 0;
                diffCoeff_next = 0;
                diffCoeff_prev = 0;
            }

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

    u = uMid;

    // // X update ------
    // interp = true;

    // // Loop over columns
    // for (int i = 1; i < mesh.nCellsX+1; i++){

    //     /* matrix components */
    //     vector<double> A_diag_Y, A_off_diag_l_Y, A_off_diag_r_Y, b_vec_Y;

    //     // Matrix construction
    //     for (int j = 1; j < mesh.nCellsY+1; j++){

    //         //obtaining cell specific coefficients
    //         diffCoeff = getCoeffFunc(uMid[i][j], i, j, sysPtr, interp);
    //         diffCoeff_next = getCoeffFunc(uMid[i+1][j], i+1, j, sysPtr, interp);
    //         diffCoeff_prev = getCoeffFunc(uMid[i-1][j], i-1, j, sysPtr, interp);
    //         // double Lr = 0.05;
    //         // double x = indexToPos(mesh, i, 'x');
    //         // double y = indexToPos(mesh, j, 'y');
    //         // double diffCoeffInterest = 1.0 * (cos(myPI * x / (2*Lr))+1) * (cos(myPI*y/(4 * Lr))+1)/4;

    //         // if (fabs(x) < Lr && fabs(y) < 2 * Lr){
    //         //     diffCoeff = diffCoeffInterest;
    //         //     diffCoeff_next = diffCoeffInterest;
    //         //     diffCoeff_prev = diffCoeffInterest;
    //         // } else {
    //         //     diffCoeff = 0;
    //         //     diffCoeff_next = 0;
    //         //     diffCoeff_prev = 0;
    //         // }

    //         // creating A matrix 
    //         A_diag_Y.push_back(1 + 2 * diffCoeff * alpha_y);
    //         A_off_diag_l_Y.push_back(-1 * diffCoeff * alpha_y);
    //         A_off_diag_r_Y.push_back(-1 * diffCoeff * alpha_y);

    //         // creating B matrix
    //         if (i == 1){
    //             if (diffBCType == "Neumann"){
    //                 b_vec_term = (1 - (2/3) * alpha_x * diffCoeff) * uMid[i][j][varIdx] +  ((2/3) * alpha_x * diffCoeff_next * uMid[i+1][j][varIdx]); 
    //             } else {
    //                 throw runtime_error("diffBCType: " + diffBCType + " is not available.");
    //             }
    //         } else if (i == mesh.nCellsX){
    //             if (diffBCType == "Neumann"){
    //                 b_vec_term = (1 - (2/3) * alpha_x * diffCoeff) * uMid[i][j][varIdx] +  ((2/3) * alpha_x * diffCoeff_prev * uMid[i-1][j][varIdx]); 
    //             } else {
    //                 throw runtime_error("diffBCType: " + diffBCType + " is not available.");
    //             }                
    //         } else {
    //             b_vec_term = (1 - 2 * alpha_x * diffCoeff) * uMid[i][j][varIdx] + alpha_x * (diffCoeff_prev * uMid[i-1][j][varIdx] + diffCoeff_next * uMid[i+1][j][varIdx]);
    //         }

    //         b_vec_Y.push_back(b_vec_term);
    //     }

    //     A_off_diag_l_Y.pop_back();
    //     A_off_diag_r_Y.pop_back(); 

    //     // Modifiying A matrix elements according to BCs
    //     if (diffBCType == "Neumann"){
    //         A_diag_Y[0] = 1 + (2.0/3.0) * alpha_y * getCoeffFunc(uMid[i][1], i, 1, sysPtr, interp);
    //         A_diag_Y.back() = 1 + (2.0/3.0)  * alpha_y * getCoeffFunc(uMid[i][mesh.nCellsY], i, mesh.nCellsY, sysPtr, interp);
    //         A_off_diag_r_Y[0] = -1 * (2.0/3.0) * alpha_y * getCoeffFunc(uMid[i][2], i, 2, sysPtr, interp);
    //         A_off_diag_l_Y.back() = -1 * (2.0/3.0) * alpha_y * getCoeffFunc(uMid[i][mesh.nCellsY-1], i, mesh.nCellsY-1, sysPtr, interp);
    //     } else {
    //         throw runtime_error("diffBCType: " + diffBCType + " is not available.");
    //     }

    //     /* solution vector, with u1 term at the zeroth index, hence need i-1 shift to match u */
    //     solution = thomasAlgo(mesh, A_diag_Y, A_off_diag_l_Y, A_off_diag_r_Y, b_vec_Y, 'y');

    //     // Applying solution to u
    //     for (int j = 1; j < mesh.nCellsY+1; j++){
    //         u[i][j][varIdx] = solution[j-1];
    //     }
    // }



    BCFunc(u, mesh, sysPtr);
    sysPtr->getEoSPtr()->cacheAll(u, mesh);
}


void CN_Diffusion_Solver(Vec2D& u, vector<int> varIdxList, GetCoeffFunc getCoeffFunc, shared_ptr<SysCalcs> sysPtr, Mesh2D& mesh, double& dt, string diffBCType, BCFunc BCFunc){

    int nx = mesh.nCellsX;
    int ny = mesh.nCellsY;

    Eigen::SparseMatrix<double> n_plus_1_matrix(nx * ny, nx * ny);
    n_plus_1_matrix.setZero();

    vector<Eigen::VectorXd> RHS_array(varIdxList.size());
    for (Eigen::VectorXd& RHS : RHS_array){
        RHS.resize(nx*ny);
        RHS.setZero();
    }

    int k; // unified/flattened index
    double r_base_x, r_base_y, rx, ry, diffCoeff;
    bool interp = true;

    r_base_x = dt / (2 * pow(mesh.dx, 2.0));
    r_base_y = dt / (2 * pow(mesh.dy, 2.0));

    // constructing the bulk of the matrix
    for (int i = 2; i < nx; i++){
        for (int j = 2; j < ny; j++){
            diffCoeff =  getCoeffFunc(u[i][j], i, j, sysPtr, interp);
            // double Lr = 0.05;
            // double x = indexToPos(mesh, i, 'x');
            // double y = indexToPos(mesh, j, 'y');
            // double diffCoeffInterest = 1.0 * (cos(myPI * x / (2*Lr))+1) * (cos(myPI*y/(4 * Lr))+1)/4;;
            // if (fabs(x) < Lr && fabs(y) < 2 * Lr){
            //     diffCoeff = diffCoeffInterest;
            // } else {
            //     diffCoeff = 0;
            // }

            rx = diffCoeff * r_base_y;
            ry = diffCoeff * r_base_y;

            // calculating unified/flattened index
            k = (i-1) * ny + j - 1; // need to subtract 1 due to 0 based indexing

            // matrix for LHS (u^n+1 matrix)
            n_plus_1_matrix.insert(k, (i-2)*ny + j-1)     = - rx;    // i-1,j
            n_plus_1_matrix.insert(k, (i-1)*ny + j-2)   = - ry;      // i, j-1   
            n_plus_1_matrix.insert(k, (i-1)*ny + j-1)     = 1 + 2*rx + 2*ry; // i,j
            n_plus_1_matrix.insert(k, (i-1)*ny + j)   = - ry;   // i,j+1
            n_plus_1_matrix.insert(k, i*ny + j-1)     = - rx;   // i+1, j


            // vector for RHS
            for (int q = 0; q < varIdxList.size(); q++){
                RHS_array[q](k) += rx * u[i-1][j][varIdxList[q]];
                RHS_array[q](k) += ry * u[i][j-1][varIdxList[q]];
                RHS_array[q](k) += (1 - 2*rx - 2*ry) * u[i][j][varIdxList[q]];
                RHS_array[q](k) += ry * u[i][j+1][varIdxList[q]];
                RHS_array[q](k) += rx * u[i+1][j][varIdxList[q]];
            }


        }
    }

    /* deal with BCs */
    if (diffBCType == "Neumann"){

        // x boundaries
        for (int j = 2; j < ny; j++){

            // left boundary (i = 1)
            diffCoeff =  getCoeffFunc(u[1][j], 1, j, sysPtr, interp);
            rx = diffCoeff * r_base_y;
            ry = diffCoeff * r_base_y;
            k = j-1;

            n_plus_1_matrix.insert(k, j-2)   = - ry;                  // i, j-1
            n_plus_1_matrix.insert(k, j-1)   = 1 - 4.0/3.0*rx + 2*rx + 2*ry; // i,j
            n_plus_1_matrix.insert(k, j)     = - ry;                  // i, j+1
            n_plus_1_matrix.insert(k, ny + j-1)   = - rx + 1.0/3.0 * rx;            // i+1, j

            for (int q = 0; q < varIdxList.size(); q++){
                RHS_array[q](k) += ry * u[1][j-1][varIdxList[q]];
                RHS_array[q](k) += (1 - 4.0/3.0*rx - 2*rx - 2*ry) * u[1][j][varIdxList[q]];
                RHS_array[q](k) += ry * u[1][j+1][varIdxList[q]];
                RHS_array[q](k) += (rx + 1.0/3.0*rx) * u[2][j][varIdxList[q]];
            }

            // right boundary (i = nx)
            diffCoeff =  getCoeffFunc(u[nx][j], nx, j, sysPtr, interp);
            rx = diffCoeff * r_base_y;
            ry = diffCoeff * r_base_y;
            k = (nx-1) * ny + j - 1;

            n_plus_1_matrix.insert(k, (nx-2)*ny + j-1)   = - rx + 1.0/3.0*rx; // i-1, j
            n_plus_1_matrix.insert(k, (nx-1)*ny + j-2)   = - ry; // i, j-1
            n_plus_1_matrix.insert(k, (nx-1)*ny + j-1)   = 1 - 4.0/3.0*rx + 2*rx + 2*ry; // i, j
            n_plus_1_matrix.insert(k, (nx-1)*ny + j)   = - ry; // i, j+1

            for (int q = 0; q < varIdxList.size(); q++){
                RHS_array[q](k) += (rx + 1.0/3.0*rx) * u[nx-1][j][varIdxList[q]];
                RHS_array[q](k) += ry * u[nx][j-1][varIdxList[q]];
                RHS_array[q](k) += (1 - 4.0/3.0*rx - 2*rx - 2*ry) * u[nx][j][varIdxList[q]];
                RHS_array[q](k) += ry * u[nx][j+1][varIdxList[q]];
            }


        }


        // y boundaries
        for (int i = 2; i < nx; i++){

            // left boundary
            int j = 1;

            diffCoeff =  getCoeffFunc(u[i][j], i, j, sysPtr, interp);
            rx = diffCoeff * r_base_y;
            ry = diffCoeff * r_base_y;
            k = (i-1) * ny + j - 1;

            n_plus_1_matrix.insert(k, (i-2)*ny + j-1)     = - rx;
            n_plus_1_matrix.insert(k, (i-1)*ny + j-1)     = 1 - 4.0/3.0*ry + 2*rx + 2*ry;
            n_plus_1_matrix.insert(k, (i-1)*ny + j)   = - ry + 1.0/3.0*ry; 
            n_plus_1_matrix.insert(k, i*ny + j-1)     = - rx;

            for (int q = 0; q < varIdxList.size(); q++){
                RHS_array[q](k) += rx * u[i-1][j][varIdxList[q]]; //i-1, j
                RHS_array[q](k) += (1 - 4.0/3.0*ry - 2*rx - 2*ry) * u[i][j][varIdxList[q]];
                RHS_array[q](k) += (ry + 1.0/3.0*ry) * u[i][j+1][varIdxList[q]];
                RHS_array[q](k) += rx * u[i+1][j][varIdxList[q]];
            }

            // right boundary
            j = ny;
            diffCoeff =  getCoeffFunc(u[i][j], i, j, sysPtr, interp);
            rx = diffCoeff * r_base_y;
            ry = diffCoeff * r_base_y;
            k = (i-1) * ny + j - 1;

            n_plus_1_matrix.insert(k, (i-2)*ny + j-1)     = - rx;    // i-1,j
            n_plus_1_matrix.insert(k, (i-1)*ny + j-2)   = - ry + 1.0/3.0*ry;      // i, j-1   
            n_plus_1_matrix.insert(k, (i-1)*ny + j-1)     = 1 - 4.0/3.0*ry + 2*rx + 2*ry; // i,j
            n_plus_1_matrix.insert(k, i*ny + j-1)     = - rx;   // i+1, j


            // vector for RHS
            for (int q = 0; q < varIdxList.size(); q++){
                RHS_array[q](k) += rx * u[i-1][j][varIdxList[q]];
                RHS_array[q](k) += (ry + 1.0/3.0*ry) * u[i][j-1][varIdxList[q]];
                RHS_array[q](k) += (1 - 4.0/3.0*ry - 2*rx - 2*ry) * u[i][j][varIdxList[q]];
                RHS_array[q](k) += rx * u[i+1][j][varIdxList[q]];
            }
 
        }

        /* handle first and last row type situations */

        // first first -----
        int i = 1;
        int j = 1;
        diffCoeff =  getCoeffFunc(u[i][j], i, j, sysPtr, interp);
        rx = diffCoeff * r_base_y;
        ry = diffCoeff * r_base_y;
        k = (i-1) * ny + j - 1;


        // matrix for LHS (u^n+1 matrix)
        n_plus_1_matrix.insert(k, (i-1)*ny + j-1)     = 1 - 4.0/3.0 * rx - 4.0/3.0*ry + 2*rx + 2*ry; // i,j
        n_plus_1_matrix.insert(k, (i-1)*ny + j)   = - ry + 1.0/3.0*ry;   // i,j+1
        n_plus_1_matrix.insert(k, i*ny + j-1)     = - rx + 1.0/3.0*rx;   // i+1, j

        // vector for RHS
        for (int q = 0; q < varIdxList.size(); q++){
            RHS_array[q](k) += (1 - 4.0/3.0*rx - 4.0/3.0*ry - 2*rx - 2*ry) * u[i][j][varIdxList[q]];
            RHS_array[q](k) += (ry + 1.0/3.0*ry) * u[i][j+1][varIdxList[q]];
            RHS_array[q](k) += (rx + 1.0/3.0*rx) * u[i+1][j][varIdxList[q]];
        }

  

        // first last ------
        i = 1; j = ny;
        diffCoeff =  getCoeffFunc(u[i][j], i, j, sysPtr, interp);
        rx = diffCoeff * r_base_y;
        ry = diffCoeff * r_base_y;
        k = (i-1) * ny + j - 1;

        n_plus_1_matrix.insert(k, (i-1)*ny + j-2)   = - ry + 1.0/3.0*rx;      // i, j-1   
        n_plus_1_matrix.insert(k, (i-1)*ny + j-1)   = 1 - 4.0/3.0 * rx - 4.0/3.0*ry + 2*rx + 2*ry; // i,j
        n_plus_1_matrix.insert(k, i*ny + j-1)     = - rx + 1.0/3.0*rx;   // i+1, j

        // vector for RHS
        for (int q = 0; q < varIdxList.size(); q++){
            RHS_array[q](k) += (ry + 1.0/3.0*ry) * u[i][j-1][varIdxList[q]];
            RHS_array[q](k) += (1 - 4.0/3.0*rx - 4.0/3.0*ry - 2*rx - 2*ry) * u[i][j][varIdxList[q]];
            RHS_array[q](k) += (rx + 1.0/3.0*rx) * u[i+1][j][varIdxList[q]];
        }


        // last first ------
        i = nx; j = 1;
        k = (i-1) * ny + j - 1;
        diffCoeff =  getCoeffFunc(u[i][j], i, j, sysPtr, interp);
        rx = diffCoeff * r_base_y;
        ry = diffCoeff * r_base_y;

        // matrix for LHS (u^n+1 matrix)
        n_plus_1_matrix.insert(k, (i-2)*ny + j-1)     = - rx + 1.0/3.0 * rx;    // i-1,j
        n_plus_1_matrix.insert(k, (i-1)*ny + j-1)     = 1 - 4.0/3.0 * rx - 4.0/3.0*ry + 2*rx + 2*ry; // i,j
        n_plus_1_matrix.insert(k, (i-1)*ny + j)   = - ry + 1.0/3.0*ry;   // i,j+1

        // vector for RHS

        for (int q = 0; q < varIdxList.size(); q++){
            RHS_array[q](k) += (rx + 1.0/3.0*rx) * u[i-1][j][varIdxList[q]];
            RHS_array[q](k) += (1 - 4.0/3.0*rx - 4.0/3.0*ry - 2*rx - 2*ry) * u[i][j][varIdxList[q]];
            RHS_array[q](k) += (ry + 1.0/3.0*ry) * u[i][j+1][varIdxList[q]];
        }



        // last row -----
        i = nx; j = ny;
        k = (i-1) * ny + j - 1;
        diffCoeff =  getCoeffFunc(u[i][j], i, j, sysPtr, interp);
        rx = diffCoeff * r_base_y;
        ry = diffCoeff * r_base_y;

        // matrix for LHS (u^n+1 matrix)
        n_plus_1_matrix.insert(k, (i-2)*ny + j-1)     = - rx + 1.0/3.0 * rx;    // i-1,j
        n_plus_1_matrix.insert(k, (i-1)*ny + j-2)   = - ry + 1.0/3.0 * ry;      // i, j-1   
        n_plus_1_matrix.insert(k, (i-1)*ny + j-1)     = 1 - 4.0/3.0 * rx - 4.0/3.0*ry + 2*rx + 2*ry; // i,j

        // vector for RHS
        for (int q = 0; q < varIdxList.size(); q++){
            RHS_array[q](k) += (rx + 1.0/3.0*rx) * u[i-1][j][varIdxList[q]];
            RHS_array[q](k) += (ry + 1.0/3.0*rx) * u[i][j-1][varIdxList[q]];
            RHS_array[q](k) += (1 - 4.0/3.0*rx - 4.0/3.0*ry - 2*rx - 2*ry) * u[i][j][varIdxList[q]];
        }


    }

    
    // cout << n_plus_1_matrix << endl;
    // cout << RHS << endl;
    // cout << "\n\n" << endl;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(n_plus_1_matrix); // Perform the factorization

    if (solver.info() != Eigen::Success) {
        // decomposition failed
        std::cerr << "Decomposition failed" << std::endl;
    }

    for (int q = 0; q < varIdxList.size(); q++){
        Eigen::VectorXd solution = solver.solve(RHS_array[q]); // Solve Ax = b

        if (solver.info() != Eigen::Success) {
            // solving failed
            std::cerr << "Solving failed for varIdx: " << varIdxList[q] << std::endl;
        }

        /* updating u with solutions */
        for (int i = 1; i < nx+1; i++){
            for (int j = 1; j < ny+1; j++){
                k = (i-1) * ny + j - 1; // need to subtract 1 due to 0 based indexing
                u[i][j][varIdxList[q]] = solution(k);
            }
        }
    }
        
    BCFunc(u, mesh, sysPtr);
    sysPtr->getEoSPtr()->cacheAll(u, mesh);

}


/**
 * Crank Nicolson solver for the diffusion equation applied to thermal conduction
*/

void CN_Conduction_Solver_ADI(Vec2D& u, vector<int> varIdxList, GetCoeffFunc getCoeffFunc, shared_ptr<SysCalcs> sysPtr, Mesh2D& mesh, double& dt, string diffBCType, BCFunc BCFunc){
    
    /* intermediate matix for after x update */
    Vec2D uMid = u;
    int varIdx = varIdxList[0];

    /* diffusion coefficient */
    double diffCoeff, diffCoeff_next, diffCoeff_prev;

    /* coefficent not containing diffCoeff*/
    double alpha_x = (dt/2.0) / (2 * pow(mesh.dx, 2.0)); 
    double alpha_y = (dt/2.0) / (2 * pow(mesh.dy, 2.0)); 

    /* whether diffusion coefficient should be reinterpreted */
    bool interp = true; 

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
            // if (i == round(mesh.nCellsX/2) && j == 1 ){
            //     cout << diffCoeff << endl;
            // }
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
                    b_vec_term = (1 - (2/3) * alpha_y * diffCoeff) * sysPtr->interp_T(u[i][j]) +  ((2/3) * alpha_y * diffCoeff_next * sysPtr->interp_T(u[i][j+1])); 
                } else {
                    throw runtime_error("diffBCType: " + diffBCType + " is not available.");
                }
            } else if (j == mesh.nCellsY){
                if (diffBCType == "Neumann"){
                    b_vec_term = (1 - (2/3) * alpha_y * diffCoeff) * sysPtr->interp_T(u[i][j]) +  ((2/3) * alpha_y * diffCoeff_prev * sysPtr->interp_T(u[i][j-1])); 
                } else {
                    throw runtime_error("diffBCType: " + diffBCType + " is not available.");
                }                
            } else {
                b_vec_term = (1 - 2 * alpha_y * diffCoeff) * sysPtr->interp_T(u[i][j]) + alpha_y * (diffCoeff_prev * sysPtr->interp_T(u[i][j-1]) + diffCoeff_next * sysPtr->interp_T(u[i][j+1]));
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

            double rho = u[i][j][0];
            double pOld = sysPtr->interp_p(u[i][j]);
            double p = pOld / sysPtr->getEoSPtr()->interp_T(rho, pOld) * solution[i-1];
            double e_new = sysPtr->getEoSPtr()->get_e(rho, p);

            u[i][j][varIdx] = sysPtr->get_KE(u[i][j]) + sysPtr->get_MagE(u[i][j]) + rho * e_new;

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
                    b_vec_term = (1 - (2/3) * alpha_x * diffCoeff) * sysPtr->interp_T(uMid[i][j]) +  ((2/3) * alpha_x * diffCoeff_next * sysPtr->interp_T(uMid[i+1][j])); 
                } else {
                    throw runtime_error("diffBCType: " + diffBCType + " is not available.");
                }
            } else if (i == mesh.nCellsX){
                if (diffBCType == "Neumann"){
                    b_vec_term = (1 - (2/3) * alpha_x * diffCoeff) * sysPtr->interp_T(uMid[i][j]) +  ((2/3) * alpha_x * diffCoeff_prev * sysPtr->interp_T(uMid[i-1][j])); 
                } else {
                    throw runtime_error("diffBCType: " + diffBCType + " is not available.");
                }                
            } else {
                b_vec_term = (1 - 2 * alpha_x * diffCoeff) * sysPtr->interp_T(uMid[i][j]) + alpha_x * (diffCoeff_prev * sysPtr->interp_T(uMid[i-1][j]) + diffCoeff_next * sysPtr->interp_T(uMid[i+1][j]));
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

            // double rho = uMid[i][j][0];
            // double pOld = sysPtr->interp_p(uMid[i][j]);
            // double p = pOld / sysPtr->getEoSPtr()->interp_T(rho, pOld) * solution[j-1];
            // double e_new = sysPtr->getEoSPtr()->get_e(rho, p);

            // u[i][j][varIdx] = sysPtr->get_KE(uMid[i][j]) + sysPtr->get_MagE(uMid[i][j]) + rho * e_new;

        }
    }

    BCFunc(u, mesh, sysPtr);
    sysPtr->getEoSPtr()->cacheAll(u, mesh);
}





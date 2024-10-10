#include "utils.hpp"


/**
 * Creates a top level vector object for 2D simulations of specified lengths
*/
Vec2D makeVec2D(int lenX, int lenY){
    // creating a Vec2D of lenX
    Vec2D vector(lenX);

    for (Vec1D& innerVec : vector){
        // resizing Vec1D to lenY
        innerVec.resize(lenY);

        for (CellVec& cellVec : innerVec){
            // resizing cellVec to cellVarNums
            cellVec.resize(cellVarsNums);
        }
    }

    return vector;
}

/**
 * Creates a top level scalar object for 2D simulations of specified lengths
*/
Scalar2D makeScalar2D(int lenX, int lenY){
    Scalar2D vector(lenX);

    for (Scalar1D& innerVec : vector){
        innerVec.resize(lenY);
    }

    return vector;
}

/**
 * set all elements to 0 in a Vec2D object
*/
void setZeros(Vec2D& u){
    for (Vec1D& Vec1D : u){
        for (CellVec& cellArray : Vec1D){
            for (double& var : cellArray){
                var = 0;
            }
        }
    }
}


/**
 * Gets a slice of a 2D matrix
*/
Scalar1D getRow(const Scalar2D& matrix, int idx){
    Scalar1D outputVec;

    if (int(matrix[1].size())-1 < idx ){
        throw invalid_argument("idx: "+to_string(idx)+" is too large for a y slice of the given matrix");
    }

    for (const Scalar1D& column : matrix){
        outputVec.push_back(column[idx]);
    }

    return outputVec;
}

/**
 * The proper way to get a columns, without swapped indices
*/
Scalar1D getColumnEoS(const Scalar2D& matrix, int idx){
    Scalar1D outputVec;

    if (int(matrix[1].size())-1 < idx ){
        throw invalid_argument("idx: "+to_string(idx)+" is too large for a y slice of the given matrix");
    }

    for (const Scalar1D& row : matrix){
        outputVec.push_back(row[idx]);
    }

    return outputVec;
}

Scalar1D getColumn(const Scalar2D& matrix, int idx){
    Scalar1D outputVec;

    if (int(matrix.size()) < idx){
        throw invalid_argument("idx: "+to_string(idx)+"is too large for a x slice of the given matrix");
    }

    outputVec = matrix[idx];
    return outputVec;
}


/**
 * Gets a slice of a 2D matrix
*/
Vec1D getRow(const Vec2D& matrix, int idx){
    Vec1D outputVec;

    if (int(matrix[1].size())-1 < idx ){
        throw invalid_argument("idx: "+to_string(idx)+" is too large for a y slice of the given matrix");
    }

    for (const Vec1D& column : matrix){
        outputVec.push_back(column[idx]);
    }

    return outputVec;
}

/* similar to range(start, end, gap) in python */
vector<int> range(int start, int end, int gap){
    vector<int> output;
    int value;
    for (int i = 0; i <= (end-start)/gap; i++){
        value = start + i * gap;
        output.push_back(value);
    }

    return output;
}

/* similar to linspace in python*/
vector<double> linspace(double start, double end, int nPoints){
    vector<double> output;
    double gap = (end - start)/(nPoints-1);

    double value;

    for (int i = 0; i < nPoints; i++){
        value = start + i * gap;
        output.push_back(value);
    }

    return output;
}

/* similar to linspace in python but samples exponentially */
vector<double> linspaceLog(double start, double end, int nPoints){
    
    double logEnd = log10(end);
    double logStart = log10(start);
    vector<double> xRange = linspace(logStart, logEnd, nPoints);
    
    vector<double> output;
    double value;

    for (int i = 0; i < nPoints; i++){
        value = pow(10, xRange[i]);
        output.push_back(value);
    }

    return output;
}

Eigen::Vector3d get_J(Vec2D& u, Mesh2D& mesh, int i, int j){
    Eigen::Vector3d J;
    double Jx, Jy, Jz;
    double dBy_dx, dBx_dy;


    // first order differences on boundaries

    /* getting Jx and dBx_dy */
    if (j == 0){
        Jx = (u[i][j+1][7] - u[i][j][7])/mesh.dy;
        dBx_dy = (u[i][j+1][5] - u[i][j][5])/mesh.dy;
    } else if (j == mesh.nCellsY+1){
        Jx = (u[i][j][7] - u[i][j-1][7])/mesh.dy;
        dBx_dy = (u[i][j][5] - u[i][j-1][5])/mesh.dy;
    } else {
        Jx = (u[i][j+1][7] - u[i][j-1][7])/(2*mesh.dy);
        dBx_dy = (u[i][j+1][5] - u[i][j-1][5])/(2*mesh.dy);
    }

    /* getting Jy and dBy_dx*/
    if (i == 0){
        Jy = -(u[i+1][j][7] - u[i][j][7])/mesh.dx;
        dBy_dx = -(u[i+1][j][6] - u[i][j][6])/mesh.dx;
    } else if (i == mesh.nCellsX+1){
        Jy = (u[i][j][7] - u[i-1][j][7])/mesh.dx;
        dBy_dx = (u[i][j][6] - u[i-1][j][6])/mesh.dx;
    } else {
        Jy = (u[i+1][j][7] - u[i-1][j][7])/(2*mesh.dx);
        dBy_dx = (u[i+1][j][6] - u[i-1][j][6])/(2*mesh.dx);
    }

    J << Jx, Jy, dBy_dx - dBx_dy;

    return J;
}

void enforce_symmetry(Vec2D& u, Mesh2D& mesh, double r_tol_x, double r_tol_y, bool enforce_x, bool enforce_y){

    if (mesh.nCellsX % 2 != 0 or mesh.nCellsY % 2 != 0){
        throw runtime_error("enforce_symmetry() only works if there are even numbers of cells in x and y");
    }
    
    /* enforcing symmetry along x direction (symmetric about the y axis) */
    double cell_left_value, cell_right_value;
    double abs_diff, percentage_diff;
    double average;
    bool is_opp_sign;

    if (enforce_x){
        for (int i = 1; i < (mesh.nCellsX/2)+1; i++){
            for (int j = 1; j < mesh.nCellsY+1; j++){
                for (int var = 0; var < cellVarsNums; var++){
                    cell_left_value = u[i][j][var];
                    cell_right_value = u[mesh.nCellsX-i][j][var];

                    if (cell_left_value != cell_right_value){
                        abs_diff = fabs(cell_left_value) - fabs(cell_right_value);
                        percentage_diff = fabs(abs_diff / cell_left_value);

                        // check if the cells have different signs
                        if (cell_left_value * cell_right_value < 0){
                            is_opp_sign = true;
                        } else {
                            is_opp_sign = false;
                        }

                        // cout << "precentage_diff: " << percentage_diff << endl;

                        if (percentage_diff != 0){
                            if (percentage_diff < r_tol_x ){
                                average = (fabs(cell_left_value) + fabs(cell_right_value))/2;
                                // average = fabs(cell_left_value);

                                if (! is_opp_sign){
                                    u[i][j][var] = average;
                                    u[mesh.nCellsX-i][j][var] = average;
                                } else {
                                    u[i][j][var] = sgn(cell_left_value) * average;
                                    u[mesh.nCellsX-i][j][var] = sgn(cell_right_value) * average;
                                }

                                // cout << "acted!" << endl;
                            }
                        }
                    
                    }
                }
        
            }
        } 
    } // end of x loop

    if (enforce_y){
        for (int j = 1; j < (mesh.nCellsY/2)+1; j++){
            for (int i = 1; i < mesh.nCellsX+1; i++){
                for (int var = 0; var < cellVarsNums; var++){
                    cell_left_value = u[i][j][var];
                    cell_right_value = u[i][mesh.nCellsY-j][var];

                    if (cell_left_value != cell_right_value){
                        abs_diff = fabs(cell_left_value) - fabs(cell_right_value);
                        percentage_diff = fabs(abs_diff / cell_left_value);

                        // check if the cells have different signs
                        if (cell_left_value * cell_right_value < 0){
                            is_opp_sign = true;
                        } else {
                            is_opp_sign = false;
                        }

                        // cout << "precentage_diff: " << percentage_diff << endl;

                        if (percentage_diff != 0){
                            if (percentage_diff < r_tol_y ){
                                average = (fabs(cell_left_value) + fabs(cell_right_value))/2;
                                // average = fabs(cell_left_value);

                                if (! is_opp_sign){
                                    u[i][j][var] = average;
                                    u[i][mesh.nCellsY-j][var] = average;
                                } else {
                                    u[i][j][var] = sgn(cell_left_value) * average;
                                    u[i][mesh.nCellsY-j][var] = sgn(cell_right_value) * average;
                                }

                                // cout << "acted!" << endl;
                            }
                        }
                    }
                    
                }
            }
        
        }
    } // end of y loop

}

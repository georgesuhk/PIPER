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


#include "mathsUtils.hpp"

// LINEAR ALGEBRA TOOLS ======

/**
 * returns the dot product between 2 vectors
*/
double dotProduct(array<double, 3> a, array<double, 3> b){
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/**
 * returns the dot product between 2 vectors
*/
double dotProduct(array<double, 2> a, array<double, 2> b){
    return a[0] * b[0] + a[1] * b[1];
}

/**
 * returns a + b
*/
array<double, 2> vectorAdd(array<double,2> a, array<double,2> b){
    array<double,2> sum;

    sum[0] = a[0] + b[0];
    sum[1] = a[1] + b[1];

    return sum; 
}

array<double, 3> vectorAdd(array<double,3> a, array<double,3> b){
    array<double,3> sum;

    sum[0] = a[0] + b[0];
    sum[1] = a[1] + b[1];
    sum[2] = a[2] + b[2];

    return sum; 
}

/**
 * returns a - b
*/
array<double, 2> vectorSubtract(array<double,2> a, array<double,2> b){
    array<double,2> sum;

    sum[0] = a[0] - b[0];
    sum[1] = a[1] - b[1];

    return sum; 
}

array<double, 3> vectorSubtract(array<double,3> a, array<double,3> b){
    array<double,3> sum;

    sum[0] = a[0] - b[0];
    sum[1] = a[1] - b[1];
    sum[2] = a[2] - b[2];

    return sum; 
}

/**
 * returns the result of a scalar multiplied with a vector
*/
array<double,2> scalarTimesVec(double scalar, array<double,2> vec){
    array<double,2> outVec;

    outVec[0] = scalar*vec[0];
    outVec[1] = scalar*vec[1];

    return outVec;
}

/**
 * returns the result of a scalar multiplied with a vector
*/
array<double,2> vecDivideByScalar(double scalar, array<double,2> vec){
    array<double,2> outVec;

    outVec[0] = vec[0] / scalar;
    outVec[1] = vec[1] / scalar;

    return outVec;
}

array<double,2> normaliseVec(array<double,2> vec){
    array<double,2> outVec;

    double vecMag = getVectorMagnitude(vec);
    outVec = vecDivideByScalar(vecMag, vec);

    return outVec;
}

/**
 * Returns the magnitude of a vector
*/
double getVectorMagnitude(array<double, 3> a){
    return sqrt(pow(a[0], 2.0) + pow(a[1], 2.0) + pow(a[2], 2.0));
}

double getVectorMagnitude(array<double, 2> a){
    return sqrt(pow(a[0], 2.0) + pow(a[1], 2.0)); 
}

/**
 * Gets the tangential unit vector to the normal vector 
*/
array<double,2> getTangVec(array<double, 2> normVec){
    double tx = normVec[1];
    double ty = -normVec[0];

    return {tx, ty};
}

/* Gets the distance to a line */
double distanceToLine(double x, double y, double diagOffset){
    // Line parameters
    double a = 1.0;
    double b = 1.0;
    double c = -diagOffset;
    
    // Calculate distance using the formula
    double numerator = a * x + b * y + c;
    double denominator = sqrt(pow(a, 2.0) + pow(b, 2.0));
    double distance = numerator / denominator;
    
    return distance;
}

/* Thomas algorithm for solving tridiagonal matrices */
vector<double> thomasAlgo(Mesh2D& mesh, vector<double> A_diag, vector<double> A_offdiag_l, vector<double> A_offdiag_r, vector<double> Bvector, char axis){

    //making all arrays the required lengths

    /* the off diagonal axis right of the main diagonal */
    A_offdiag_r.push_back(0.0);

    /* the off diagonal axis left of the main diagonal */
    A_offdiag_l.insert(A_offdiag_l.begin(), 0.0);

    vector<double> fArray = {0};
    vector<double> deltaArray = {0};

    double nCells;
    if (axis == 'x'){
        nCells = mesh.nCellsX;
    } else if (axis == 'y'){
        nCells = mesh.nCellsY;
    } else {
        throw runtime_error("axis invalid in thomasAlgo()");
    }

    // forward elimination
    double f, delta;
    for (int i = 0; i < nCells; i++){
        f = A_offdiag_r[i] / (A_diag[i] - A_offdiag_l[i] * fArray[i]);
        delta = (Bvector[i] - A_offdiag_l[i] * deltaArray[i]) / (A_diag[i] - A_offdiag_l[i] * fArray[i]);

        fArray.push_back(f);
        deltaArray.push_back(delta);
    }

    fArray.push_back(0); // last element of fArray should always be 0
    vector<double> uArray(nCells+1);

    // back subtitution
    for (int i = nCells-1; i >= 0; --i){
        uArray[i] = deltaArray[i+1] - fArray[i+1] * uArray[i+1];
    }

    vector<double> solution(uArray.begin(), uArray.end()-1);
    return solution;
}






// OTHER MATHS TOOLS ======

// simple functions ------
/* sgn function, returns the sign of the input */
int sgn(double x){
    if (x >= 0){
        return 1;
    } else {
        return -1;
    }
}

int sgn(int x){
    if (x >= 0){
        return 1;
    } else {
        return -1;
    }
}

double minBetween(double a, double b){
    double output;
    if (a < b){
        output = a; 
    } else {
        output = b;
    }
    return output;
}

int minBetween(int a, int b){
    int output;
    if (a < b){
        output = a; 
    } else {
        output = b;
    }
    return output;
}

double maxBetween(double a, double b){
    double output;
    if (a > b){
        output = a; 
    } else {
        output = b;
    }
    return output;
}

int maxBetween(int a, int b){
    int output;
    if (a > b){
        output = a; 
    } else {
        output = b;
    }
    return output;
}

/* return equal if two numbers are within a tolerance */
bool floatsAreClose(double a, double b, double tol){
    return fabs(a - b) <= tol;
}

double roundToPlace(double value, double roundFactor){
    return round(value * roundFactor)/roundFactor;
}



// advanced functions ------


/**
 * Binlinear interpolaton of cellArray values to obtain cellArray at given xPos and yPos
*/
CellVec bilinearInterp(const Vec2D& u, Mesh2D& mesh, double xPos, double yPos){
    CellVec interpResults;


    int i = round(((xPos / mesh.dx) - mesh.xMin) + 1);
    int j = round(((yPos / mesh.dy) - mesh.yMin) + 1);

    // cout << "i: " << i << " j: " << j << endl;


    double x1 = mesh.xMin + (i-2) * mesh.dx;
    double x2 = mesh.xMin + (i) * mesh.dx;
    double y1 = mesh.yMin + (j-2) * mesh.dy;
    double y2 = mesh.yMin + (j) * mesh.dy;

    //coefficients used within the interpolation in each direction
    double xInterpCoeff1 = (x2 - xPos) / (x2 - x1);
    double xInterpCoeff2 = (xPos - x1) / (x2 - x1);

    double yInterpCoeff1 = (y2 - yPos) / (y2 - y1);
    double yInterpCoeff2 = (yPos - y1) / (y2 - y1);

    //for each cellArray var
    double xInterp1, xInterp2, interpResult;

    for (int var = 0; var < (cellVarsNums-1); var++){
        if (var == 0){
            // cout << "u_x1 = " << u[i-1][j-1][var] << ", u_x2: " << u[i+1][j-1][var] << endl;
        }
        xInterp1 = xInterpCoeff1 * u[i-1][j-1][var] + xInterpCoeff2 * u[i+1][j-1][var];
        xInterp2 = xInterpCoeff1 * u[i-1][j+1][var] + xInterpCoeff2 * u[i+1][j+1][var];

        interpResults[var] = yInterpCoeff1 * xInterp1 + yInterpCoeff2 * xInterp2;
        //override
        interpResults[var] = u[i][j][var];
    }

    // cout << "rho interp: " << interpResults[0] << endl;
    // cout << "rho no interp: " << u[i][j][0] << endl;

    return interpResults;
}


double bilinearInterp(vector<vector<double>>& data, double& xVar, double& yVar, double& x1, double& x2, double xLowerIdx,
double& y1, double& y2, double yLowerIdx, bool verbose){
    double interpResults;

    //coefficients used within the interpolation in each direction
    double x1_coeff = (x2 - xVar) / (x2 - x1);
    double x2_ceoff = (xVar - x1) / (x2 - x1);

    double y1_coeff = (y2 - yVar) / (y2 - y1);
    double y2_coeff = (yVar - y1) / (y2 - y1);

    double q_11 = data[yLowerIdx][xLowerIdx];
    double q_21 = data[yLowerIdx][xLowerIdx+1];
    double q_12 = data[yLowerIdx+1][xLowerIdx];
    double q_22 = data[yLowerIdx+1][xLowerIdx+1];

    //for each cellArray var
    double q_p1, q_p2;

    q_p1 = x1_coeff * q_11 + x2_ceoff * q_21;
    q_p2 = x1_coeff * q_12 + x2_ceoff * q_22;

    if (verbose){
        cout << "y1_coeff: " << y1_coeff << endl;
        cout << "y2_coeff: " << y2_coeff << endl;

        cout << "q_p1: " << q_p1 << endl;
        cout << "q_p2: " << q_p2 << endl;
    }




    interpResults = y1_coeff * q_p1 + y2_coeff * q_p2;

    return interpResults;
}
/**
 * uses the Bisect algorithm to find the lower index that bounds the value searched from in a vector of data
 * activeRange specifies the indices to set as initial Bisection bounds
*/
int getLowerBound(const double& val, array<int,2> activeRange, vector<double>& data){
    int idxNeg = activeRange[0];
    int idxPos = activeRange[1];
    int idxMid;

    double valMid, fMid;

    int steps = 0;
    while (fabs(idxPos - idxNeg) > 1){
        steps += 1;
        idxMid = round((idxPos + idxNeg)/2);

        valMid = data[idxMid];
        fMid = valMid - val;

        if (fMid < 0){
            idxNeg = idxMid;
        } else {
            idxPos = idxMid;
        }   
    } 

    // cout << "steps taken to find val: " << steps << endl;
    return idxNeg;
}


/** 
 * Bisection Root Finder used in TabEoS
*/
double BisectSolver(double& rho, int& rhoLowerIdx, double& e, Scalar1D& rhoData, Scalar2D& eData, Scalar1D pData, double atol, int maxSteps, bool verbose){
    double p;

    int rhoHigherIdx = rhoLowerIdx+1;
    double rhoLower = rhoData[rhoLowerIdx];
    double rhoHigher = rhoData[rhoHigherIdx];

    // coefficients for bilinear interpolation
    double rho1_coeff = (rhoHigher - rho) / (rhoHigher - rhoLower);
    double rho2_coeff = (rho - rhoLower) / (rhoHigher - rhoLower);

    double rho_lin_coeff = (rho - rhoLower)/(rhoHigher - rhoLower);


    double idxNeg = 0;
    double idxPos = pData.size()-1;
    int idxMid;

    double eMid, fMid;
    int step = 0;

    while (fabs(idxPos - idxNeg) > 1 && step < maxSteps){
        step += 1;
        //compute midpoint
        idxMid = round((idxPos + idxNeg)/2);

        //more accurate way
        eMid = eData[idxMid][rhoLowerIdx] + (eData[idxMid][rhoHigherIdx] - eData[idxMid][rhoLowerIdx]) * rho_lin_coeff;
        fMid = eMid - e;

        if (fMid < 0){
            idxNeg = idxMid;
        } else {
            idxPos = idxMid;
        }
    }

    // idxNeg = idxNeg;
    // idxPos = idxPos;

    step = 0;
    double pGap = pData[idxPos] - pData[idxNeg];
    double pLower = pData[idxNeg];
    double pHigher = pData[idxPos];
    double pMid;

    // subscripts here are for rho, then p (opposite to the actual indexing)
    double e_11 = eData[idxNeg][rhoLowerIdx];
    double e_21 = eData[idxNeg][rhoHigherIdx];
    double e_12 = eData[idxPos][rhoLowerIdx];
    double e_22 = eData[idxPos][rhoHigherIdx];

    double eMid_p1 = rho1_coeff * e_11 + rho2_coeff * e_21;
    double eMid_p2 = rho1_coeff * e_12 + rho2_coeff * e_22;

    double p1_coeff, p2_coeff;

    // convergence criterion: can use (e - eMid) or (pHigher - pLower)
    while (fabs(pHigher - pLower) < atol && step < maxSteps){
        step += 1;
        pMid = (pLower + pHigher)/2;
        
        // cout << "step: " << step << endl;
        
        // y linear interpolate coeffs
        p1_coeff = (pData[idxPos] - pMid) / pGap;
        p2_coeff = (pMid - pData[idxNeg]) / pGap;

        eMid = p1_coeff * eMid_p1 + p2_coeff * eMid_p2;

        fMid = eMid - e;
        // cout << "ediff: " << eMid - e << endl;

        if (fMid < 0){
            pLower = pMid;
        } else {
            pHigher = pMid;
        }
    }

    p = (pLower + pHigher)/2;

    if (verbose){
        cout << "e: " << e << endl;
        cout << "eDiff: " << eMid - e << endl;
        cout << "p interp lower idx: " << idxNeg << endl;
        cout << "p1_coeff: " << p1_coeff << endl;
        cout << "p2_coeff: " << p2_coeff << endl;
        cout << "eMid_p1: " << eMid_p1 << endl;
        cout << "eMid_p2: " << eMid_p2 << endl;
    }




    return p;
}

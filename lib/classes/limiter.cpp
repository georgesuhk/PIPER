#include "limiter.hpp"

/* outer function for limiters that returns a top level vector with limiter values for each cell */
Vec2D getLimiterVec(const Vec2D& u, Mesh2D& mesh, Limiter limiter, char axis){
    Vec2D limiterVec = makeVec2D(mesh.nCellsX+1, mesh.nCellsY+1); //ghost cell at start that shold never be accessed
    double r;

    for (int i = 1; i < mesh.nCellsX+1; i++){
        for (int j = 1; j < mesh.nCellsY+1; j++){
            for (int var = 0; var < (cellVarsNums); var++){
                //slope assessment variable
                if (axis == 'x'){
                    r = (u[i][j][var] - u[i-1][j][var])/(u[i+1][j][var] - u[i][j][var]);
                    if ( (u[i+1][j][var] - u[i][j][var]) == 0 ){
                        r = numeric_limits<double>::max();
                    }
                } else if (axis == 'y'){
                    r = (u[i][j][var] - u[i][j-1][var])/(u[i][j+1][var] - u[i][j][var]);
                    if ( (u[i][j+1][var] - u[i][j][var]) == 0 ){
                        r = numeric_limits<double>::max();
                    }                       
                } else {
                    throw invalid_argument("In getLimiterVec(): axis argument invalid");
                }

                limiterVec[i][j][var] = limiter(r);
            }
        }
    }

    return limiterVec;
}



/**
 * Minbee limiter function
 * r is the slope assessment variable
*/
double minbee(double r){
    double output;
    double eR;

    if (r <= 0){
        output = 0;
    } 
    else if (r <= 1){
        output = r;
    } else {
        eR = 2 / (1 + r);
        output = minBetween(1.0, eR);
    }
    return output;
}

/**
 * Van leer limiter function
*/
double vanLeer(double r){
    double output;
    double eR;

    if (r <= 0){
        output = 0;
    } 
    else {
        eR = 2 / (1 + r);
        output = minBetween( (2*r)/(1+r) , eR);        
    }
    return output;
}



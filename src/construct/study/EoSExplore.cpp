#include "mutation++.h"
#include "simulation.hpp"

using namespace Mutation;
using namespace Mutation::Thermodynamics;
using namespace Mutation::Transport; 
using namespace std;
 
string mixName = "HFusion";
// double rho = 0.01;

double p = 101325;
double temp = 100000;
const int cellVarsNums = 9;

int main(){

    Mixture mix(mixName);
    mix.equilibrate(temp, p);

    cout << "num density: " << mix.numberDensity() << endl;


    cout << "no. of collision pairs: " << mix.nCollisionPairs() << endl;    

    CollisionDB ColDB = mix.collisionDB();

    cout << "no. of collision species " << ColDB.nSpecies() << endl;

    cout << "Q11 (e - e collision integral): " << ColDB.Q11ee() << endl;

    Eigen::ArrayXd Q11ei = ColDB.Q11ei();

    cout << "Q11 (e - i collisions): " << endl;
    cout << Q11ei << endl;

    Eigen::ArrayXd Q11ii = ColDB.Q11ii();    

    cout << "Q11 (i - i collisions): " << endl;
    cout << Q11ii << endl;

    // need to find meaning of Q11, Q12, Q13...


    // diffusion coefficients

    Eigen::ArrayXd ei_diff_coeff = ColDB.nDei();
    cout << "diff coefficient ei: " << endl;
    cout << ei_diff_coeff << endl;

    return 0; 
}
#ifndef EOSUTILS_HPP
#define EOSUTILS_HPP

#include "calcs.hpp"
#include "utils.hpp"
#include "mutation++.h"

using namespace Mutation;
using namespace Mutation::Thermodynamics;
using namespace Mutation::Transport; 

double BisectSolverEoS_T(double& rho, double& p, double& particleMass, Mixture& mix, double atol, double rtol, int maxSteps, double searchTol);

double BisectSolverEoS_P(double& rho, double& temp, double& particleMass, Mixture& mix, double atol, double rtol, int maxSteps, double searchTol);

#endif

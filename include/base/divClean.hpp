#ifndef DIVCLEAN_HPP
#define DIVCLEAN_HPP

#include "utils.hpp"
#include "sysCalcs.hpp"
#include "BC.hpp"

/* Getting parabolic wave speed based on an empirically determined factor */
double getCp(double ch);

/* The parabolic source term update for psi */
double paraUpdatePsi(double psiOld, double ch, double dt);

/* return the flux for the separated GLM divergence cleaning system */
array<double,2> getDCFlux(const CellVec& uL, const CellVec& uR, double ch, char axis);

#endif
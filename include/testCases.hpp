#ifndef TESTCASES_HPP
#define TESTCASES_HPP

#include "settings.hpp"

extern double rho_SF;
extern double p_SF;
double wx = 0.0;

//BW Test

vector<CellVec> BrioWuTestX = {{1, 0, 0, 0, 1, 0.75, 1, 0, 0}, {0.125, 0, 0, 0, 0.1, 0.75, -1, 0, 0}};
vector<CellVec> BrioWuTestX_SI = {{1*rho_SF, 0, 0, 0, 1*p_SF, 0.75*sqrt(p_SF), 1*sqrt(p_SF), 0, 0}, {0.125*rho_SF, 0, 0, 0, 0.1*p_SF, 0.75*sqrt(p_SF), -1*sqrt(p_SF), 0, 0}};

vector<CellVec> BrioWuPIP = {{1*rho_SF, 0, 0, 0, 1*p_SF, 0.75*sqrt(p_SF), 1*sqrt(p_SF), 0, 0, wx*sqrt(p_SF), 0, 0}, {0.125*rho_SF, 0, 0, 0, 0.1*p_SF, 0.75*sqrt(p_SF), -1*sqrt(p_SF), 0, 0, wx*sqrt(p_SF), 0, 0}};

vector<CellVec> Toro1PIP = {{1*rho_SF, 0, 0, 0, 1*p_SF, 0, 0, 0, 0, wx*sqrt(p_SF), 0, 0}, {0.125*rho_SF, 0, 0, 0, 0.1*p_SF, 0, 0, 0, 0, wx*sqrt(p_SF), 0, 0}};

// vector<CellVec> BrioWuTestX_SI = {{1*rho_SF, 0, 0, 0, 1*p_SF*pAtmos, 0.75*sqrt(pAtmos*p_SF), 1*sqrt(pAtmos*p_SF), 0, 0}, {0.125*rho_SF, 0, 0, 0, 0.1*p_SF*pAtmos, 0.75*sqrt(pAtmos*p_SF), -1*sqrt(pAtmos*p_SF), 0, 0}};
// vector<cellArray> BrioWuTestX_SINoScale = {{1*rho_SF, 0, 0, 0, 1*p_SF*pAtmos, 0.75*sqrt(pAtmos*p_SF)*sqrt(vacPermeab), 1*sqrt(pAtmos*p_SF)*sqrt(vacPermeab), 0, 0}, {0.125*rho_SF, 0, 0, 0, 0.1*p_SF*pAtmos, 0.75*sqrt(pAtmos*p_SF)*sqrt(vacPermeab), -1*sqrt(pAtmos*p_SF)*sqrt(vacPermeab), 0, 0}};


vector<CellVec> BW_NoB_SI = {{1*rho_SF, 0, 0, 0, 1*p_SF*pAtmos, 0, 0, 0, 0}, {0.125*rho_SF, 0, 0, 0, 0.1*p_SF*pAtmos, 0, 0, 0, 0}};





#endif




#ifndef TESTCASES_HPP
#define TESTCASES_HPP

#include "settings.hpp"

extern double rho_SF;
extern double p_SF;
double wx = 0.0;

//BW Test

vector<CellVec> BrioWuTestX = {{1, 0, 0, 0, 1, 0.75, 1, 0, 0}, {0.125, 0, 0, 0, 0.1, 0.75, -1, 0, 0}};
vector<CellVec> BrioWuTestX_SI = {{1*rho_SF, 0, 0, 0, 1*p_SF, 0.75*sqrt(p_SF), 1*sqrt(p_SF), 0, 0}, {0.125*rho_SF, 0, 0, 0, 0.1*p_SF, 0.75*sqrt(p_SF), -1*sqrt(p_SF), 0, 0}};


// PIP shock tube tests
vector<CellVec> BrioWuPIP = {{1*rho_SF, 0, 0, 0, 1*p_SF, 0.75*sqrt(p_SF), 1*sqrt(p_SF), 0, 0, wx*sqrt(p_SF), 0, 0}, {0.125*rho_SF, 0, 0, 0, 0.1*p_SF, 0.75*sqrt(p_SF), -1*sqrt(p_SF), 0, 0, wx*sqrt(p_SF), 0, 0}};

// ran to 0.2
double RJ_wx = -3.0;
vector<CellVec> RyuJonesPIP = {{1.08*rho_SF, 1.2*sqrt(p_SF)/sqrt(rho_SF), 0.01*sqrt(p_SF)/sqrt(rho_SF), 0.5*sqrt(p_SF)/sqrt(rho_SF), 0.95*p_SF, 0.564190*sqrt(p_SF), 1.015541*sqrt(p_SF), 0.564190*sqrt(p_SF), 0, RJ_wx*sqrt(p_SF)/sqrt(rho_SF), 0, 0}, {0.9891*rho_SF, -0.0131*sqrt(p_SF)/sqrt(rho_SF), 0.0269*sqrt(p_SF)/sqrt(rho_SF), 0.010037*sqrt(p_SF)/sqrt(rho_SF), 0.97159*p_SF, 0.564190*sqrt(p_SF), 1.135262*sqrt(p_SF), 0.564923*sqrt(p_SF), 0, RJ_wx*sqrt(p_SF)/sqrt(rho_SF), 0, 0}};
vector<CellVec> RP4PIP = {{1.0*rho_SF, 0, 0, 0, 1.0*p_SF, 1.3*sqrt(p_SF), 1.0*sqrt(p_SF), 0*sqrt(p_SF), 0, 0, 0, 0}, {0.4*rho_SF, 0, 0, 0, 0.4*p_SF, 1.3*sqrt(p_SF), -1.0*sqrt(p_SF), 0*sqrt(p_SF), 0, 0, 0, 0}};

// 2 fluid tests
double ss_B0 = 1;
double ss_P0 = 0.15;
vector<CellVec> slow_shock_PIP = {{1.0*rho_SF, 0, 0, 0, ss_P0*p_SF, 0.3*ss_B0*sqrt(p_SF), ss_B0*sqrt(p_SF), 0*sqrt(p_SF), 0, 0, 0, 0}, {1.0*rho_SF, 0, 0, 0, ss_P0*p_SF, 0.3*ss_B0*sqrt(p_SF), -ss_B0*sqrt(p_SF), 0*sqrt(p_SF), 0, 0, 0, 0}};



vector<CellVec> Toro1PIP = {{1*rho_SF, 0, 0, 0, 1*p_SF, 0.75*sqrt(p_SF), 0, 0, 0, wx*sqrt(p_SF), 0, 0}, {0.125*rho_SF, 0, 0, 0, 0.1*p_SF, 0.75*sqrt(p_SF), 0, 0, 0, wx*sqrt(p_SF), 0, 0}};

// vector<CellVec> BrioWuTestX_SI = {{1*rho_SF, 0, 0, 0, 1*p_SF*pAtmos, 0.75*sqrt(pAtmos*p_SF), 1*sqrt(pAtmos*p_SF), 0, 0}, {0.125*rho_SF, 0, 0, 0, 0.1*p_SF*pAtmos, 0.75*sqrt(pAtmos*p_SF), -1*sqrt(pAtmos*p_SF), 0, 0}};
// vector<cellArray> BrioWuTestX_SINoScale = {{1*rho_SF, 0, 0, 0, 1*p_SF*pAtmos, 0.75*sqrt(pAtmos*p_SF)*sqrt(vacPermeab), 1*sqrt(pAtmos*p_SF)*sqrt(vacPermeab), 0, 0}, {0.125*rho_SF, 0, 0, 0, 0.1*p_SF*pAtmos, 0.75*sqrt(pAtmos*p_SF)*sqrt(vacPermeab), -1*sqrt(pAtmos*p_SF)*sqrt(vacPermeab), 0, 0}};


vector<CellVec> BW_NoB_SI = {{1*rho_SF, 0, 0, 0, 1*p_SF*pAtmos, 0, 0, 0, 0}, {0.125*rho_SF, 0, 0, 0, 0.1*p_SF*pAtmos, 0, 0, 0, 0}};





#endif




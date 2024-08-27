#include "init.hpp"
#include "testCases.hpp"
#include "simulation.hpp"

// data folders ------

string dataFolder = "./EoSData/HFusion4_S1/";
string resultsFolder = "./output/MHDEoS/";

// choosing system ------

const int cellVarsNums = 12;
const int omp_threads = 1;
// const int cellVarsNums = 9;

double mass = protonMass;

// setting grid ------

int nCellsX = 350, nCellsY = 2;
double xMin = 0, xMax = 1.0;
double yMin = 0, yMax = 0.1;

double tMin = 0, tMax = 0.1*sqrt(rho_SF)/sqrt(p_SF);
int maxSteps = 1e6;
Mesh2D mesh(xMin, xMax, nCellsX, yMin, yMax, nCellsY);

// initial conditions ------

double rho_SF = 1e-9;
double p_SF = 1e-4;
vector<double> interfacePositions = {0.5};
vector<CellVec> initCellVecs = BrioWuPIP;

// source terms ------
ExplicitSolver explicitSolver = RK4;
bool doSourceUpdate = true;
int sourceTimeRatio = 1;
int impExRatio = 1;
vector<implicitSource> implicitSources = {conduction};
// vector<implicitSource> implicitSources = {};
vector<SourceFuncEx> exSourceFuncs = {};
// vector<SourceFuncEx> exSourceFuncs = {w_evolution_func};

// EoS ------
double gammaFac = 5.0/3.0;
double constResis = 3;
double constThermCon = 0;
double mass_frac_n = 0.8;
double mass_frac_i = 1.0 - mass_frac_n;

// divergence cleaning ------
bool doDC = true;

// evolver & BCs ------

BCFunc BC = TransBCs;
SLICEvolver evolver(BC, mesh);
shared_ptr<Evolver> evolverPtr = make_shared<SLICEvolver>(evolver);

// recorder and exporter ------
double recordingDelayTime = (tMax-tMin)/20;
Exporter exporter = ExportPIP;
double simExportDelay = (tMax - tMin)/20;
bool doInSimExport = true;

// forcing initial time steps ------

/* the step number up until which time step is forced */
int forced_step_lim = 0;

/* the ratio to lower time step by */
double forced_ratio = 1e-5;

int main(void){
    // SET UP ======
    omp_set_dynamic(0); 

    // IdealEoS EoSIdeal(gammaFac, mass);
    // EoSIdeal.set_constResis(constResis);
    // EoSIdeal.set_constThermCon(constThermCon);
    // EoSIdeal.set_mass_frac_n(mass_frac_n);
    // EoSIdeal.set_mass_frac_i(mass_frac_i);
    // shared_ptr<EoS> EoSPtr = make_shared<IdealEoS>(EoSIdeal);

    TabEoS EoSTab;
    EoSTab.genFromData(mesh, {"pressure","densities","T","Cs","e","gamma","resis","thermCon","mass_frac_e","mass_frac_n","mass_frac_i"}, dataFolder, ',');
    shared_ptr<EoS> EoSPtr = make_shared<TabEoS>(EoSTab);

    PIP0_Calcs sysCalcs(EoSPtr);
    shared_ptr<SysCalcs> sysPtr = make_shared<PIP0_Calcs>(sysCalcs);

    // FIPCalcs sysCalcs(EoSPtr);
    // shared_ptr<SysCalcs> sysPtr = make_shared<FIPCalcs>(sysCalcs);

    Vec2D uInit = initPlanar(initCellVecs, interfacePositions, mesh, sysPtr, BC, 'x');
    Recorder recorder;
    recorder.setDelayTime(recordingDelayTime);
    shared_ptr<Recorder> recorderPtr = make_shared<Recorder>(recorder);

    Simulation sim(uInit, evolverPtr, sysPtr, recorderPtr, BC, exporter, mesh, resultsFolder, tMax);
    sim.setDoSourceUpdate(doSourceUpdate);
    sim.setSourceTimeRatio(sourceTimeRatio, impExRatio);
    sim.setExplicitSolver(explicitSolver);
    sim.setImplicitSources(implicitSources);
    sim.setExplicitSourceFuncs(exSourceFuncs);
    sim.setDoDC(doDC);

    sim.forceRecordAll();
    sim.enableProgressUpdate(0.1);
    sim.setDoInSimExport(doInSimExport);
    sim.setExportGap(simExportDelay);
    sim.inform();

    // testing

    // double rho = 1.67e-8;
    // double e = 26900;

    // double p = sysPtr->getEoSPtr()->interp_p(rho, e);
    // double e_interp = sysPtr->getEoSPtr()->get_e(rho, p);

    // cout << "p: " << p << endl;
    // cout << "e_interp: " << e_interp << endl;


    // // SIMULATION ======

    // cout << "Starting simulation. Tmax = " << tMax << "\n\n" << endl;
    // double t = tMin, step = 1;

    // while (t <= tMax && step <= maxSteps){
        
    //     // forcing smaller time steps in the beginning
    //     if (step < forced_step_lim){
    //         double forced_dt = forced_ratio * sim.get_min_dt();
    //         sim.force_set_dt(forced_dt);
    //     } else {
    //         sim.update_dt();
    //     }
    //     sim.evolve();

    //     t = sim.getTime();
    //     step = sim.getStep();
    //     cout << "t: " << t << endl;
    //     // cout << "dt: " << sim.get_dt() << endl;
    // }

    // cout << "Simulation Completed." << endl;
    // sim.forceRecordAll();
    // sim.exportAll(resultsFolder);

    return 0; 
}
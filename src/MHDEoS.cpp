#include "init.hpp"
#include "testCases.hpp"
#include "simulation.hpp"

// data folders ------

string dataFolder = "./EoSData/HFusion4_1/";
string resultsFolder = "./output/MHDEoS/";

// choosing system ------

const int cellVarsNums = 9;
double gammaFac = 5.0/3.0;
double mass = 1.67e-27;

// setting grid ------

int nCellsX = 300, nCellsY = 4;
double xMin = 0, xMax = 1;
double yMin = 0, yMax = 0.01;

double tMin = 0, tMax = 0.1*sqrt(rho_SF)/sqrt(p_SF*pAtmos);
int maxSteps = 0;

// initial conditions ------

double rho_SF = 1.0;
double p_SF = 1.0;
vector<double> interfacePositions = {0.5};
vector<CellVec> initCellVecs = BrioWuTestX;

// source terms ------
bool doSourceUpdate = true;
int sourceTimeRatio = 1;
vector<implicitSource> implicitSources = {ohmic_diffusion};
double constResis = 1000.0;

// divergence cleaning ------
bool doDC = true;

// evolver & BCs ------

BCFunc BC = TransBCs;
SLICEvolver evolver(BC);
shared_ptr<Evolver> evolverPtr = make_shared<SLICEvolver>(evolver);

// recorder and exporter ------
double recordingDelayTime = (tMax-tMin)/0.1;
Exporter exporter = ExportFIP;

int main(void){
    // SET UP ======

    Mesh2D mesh(xMin, xMax, nCellsX, yMin, yMax, nCellsY);

    // IdealEoS EoSIdeal(gammaFac, mass);
    // EoSIdeal.set_constResis(constResis);
    // shared_ptr<EoS> EoSPtr = make_shared<IdealEoS>(EoSIdeal);

    TabEoS EoSTab;
    EoSTab.genFromData(mesh, {"pressure","densities","T","Cs","e","resis","thermCon","n_e","n_n","n_i"}, dataFolder, ',');
    shared_ptr<EoS> EoSPtr = make_shared<TabEoS>(EoSTab);

    FIPCalcs sysCalcs(EoSPtr);
    shared_ptr<SysCalcs> sysPtr = make_shared<FIPCalcs>(sysCalcs);

    Vec2D uInit = initPlanar(initCellVecs, interfacePositions, mesh, sysPtr, BC, 'x');
    Recorder recorder;
    recorder.setDelayTime(recordingDelayTime);
    shared_ptr<Recorder> recorderPtr = make_shared<Recorder>(recorder);

    Simulation sim(uInit, evolverPtr, sysPtr, recorderPtr, BC, exporter, mesh, resultsFolder, tMax);
    sim.setDoSourceUpdate(doSourceUpdate);
    sim.setSourceTimeRatio(sourceTimeRatio);
    sim.setImplicitSources(implicitSources);
    sim.setDoDC(doDC);

    sim.forceRecordAll();
    sim.enableProgressUpdate(0.05);
    sim.inform();

    // SIMULATION ======

    cout << "Starting simulation. \n\n" << endl;
    double t = tMin, step = 1;

    while (t <= tMax && step <= maxSteps){
        sim.evolve();

        t = sim.getTime();
        step = sim.getStep();
        // cout << "t: " << t << endl;
    }

    cout << "Simulation Completed." << endl;
    sim.forceRecordAll();
    sim.exportAll(resultsFolder);

    return 0; 
}
#include "init.hpp"
#include "testCases.hpp"
#include "simulation.hpp"

// data folders ------

string dataFolder = "./EoSData/HFusion4_S1/";
string rootResultsFolder = "./output/MHDEoS/stored/report/RJ_drift/";
vector<string> resultsFolderArray = {
    "cfl0.4_400cells_rho1.5e-7/",
    "cfl0.4_400cells_rho1.6-e7/",
    "cfl0.4_400cells_rho1.7e-7/",
    "cfl0.4_400cells_rho1.8e-7/",
    "cfl0.4_400cells_rho1.85e-7/"
};

// choosing system ------

const int cellVarsNums = 12;
const int omp_threads = 1;

double mass = 1.67e-27;

// setting grid ------


int nCellsX = 10, nCellsY = 2;
double xMin = 0, xMax = 1.0;
double yMin = 0, yMax = 0.1;

double tMin = 0;
int maxSteps = 1e6;
Mesh2D mesh(xMin, xMax, nCellsX, yMin, yMax, nCellsY);

// initial conditions ------

vector<double> rho_SF_Array = {1.5e-7, 1.6e-7, 1.7e-7, 1.8e-7, 1.85e-7};
int numSims = 5;
double rho_SF = 1.6e-7;
double p_SF = 1e-4;
vector<double> interfacePositions = {0.5};

// source terms ------
ExplicitSolver explicitSolver = RK4;
bool doSourceUpdate = true;
int sourceTimeRatio = 1;
int impExRatio = 5;
// vector<implicitSource> implicitSources = {ohmic_diffusion};
vector<implicitSource> implicitSources = {};
// vector<SourceFuncEx> exSourceFuncs = {};
vector<SourceFuncEx> exSourceFuncs = {w_evolution_func};

// EoS ------
double gammaFac = 5.0/3.0;
double constResis = 3;
double mass_frac_n = 0.8;
double mass_frac_i = 1.0 - mass_frac_n;

// divergence cleaning ------
bool doDC = true;

// evolver & BCs ------

BCFunc BC = TransBCs;
SLICEvolver evolver(BC, mesh);
shared_ptr<Evolver> evolverPtr = make_shared<SLICEvolver>(evolver);

// recorder and exporter ------
Exporter exporter = ExportPIP;
bool doInSimExport = false;

// forcing initial time steps ------

/* the step number up until which time step is forced */
int forced_step_lim = 10;

/* the ratio to lower time step by */
double forced_ratio = 1e-3;

int main(void){
    // SET UP ======
    omp_set_dynamic(0); 

    // IdealEoS EoSIdeal(gammaFac, mass);
    // EoSIdeal.set_constResis(constResis);
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

    for (int simNum = 0; simNum < numSims; simNum++){
        double rho_SF = rho_SF_Array[simNum];
        string resultsFolder = rootResultsFolder+resultsFolderArray[simNum];

        double tMax = 0.2*sqrt(rho_SF)/sqrt(p_SF);
        double recordingDelayTime = (tMax-tMin)/20;
        double simExportDelay = (tMax - tMin)/20;

        vector<CellVec> initCellVecs = RyuJonesPIP;

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

        // SIMULATION ======

        cout << "Starting simulation. Tmax = " << tMax << "\n\n" << endl;
        double t = tMin, step = 1;

        while (t <= tMax && step <= maxSteps){
            
            // forcing smaller time steps in the beginning
            if (step < forced_step_lim){
                double forced_dt = forced_ratio * sim.get_min_dt();
                sim.force_set_dt(forced_dt);
            } else {
                sim.update_dt();
            }
            sim.evolve();

            t = sim.getTime();
            step = sim.getStep();
        }

        cout << "Simulation Completed." << endl;
        sim.forceRecordAll();
        sim.exportAll(resultsFolder);
     
    }

    return 0; 
}
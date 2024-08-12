#include "simulation.hpp"
#include "detachInit.hpp"

using namespace std;
const int cellVarsNums = 12;

// grid ------

int nCellsX = 30, nCellsY = 2;
double xMin = 0, xMax = 1;
double yMin = 0, yMax = 0.2;
// double tMin = 0, tMax = 0.1/sqrt(pAtmos);
double tMin = 0, tMax = 15e-7;

// double tMin = 0, tMax = 200e-7;


Mesh2D mesh(xMin, xMax, nCellsX, yMin, yMax, nCellsY);


// init params ------

const int omp_threads = 1;
double rho_SF = 1.8e-7;
double p_SF = 1e-4;
vector<double> initParams = {1*p_SF, 200000, 5000, 1*sqrt(p_SF), 0};
int maxSteps = 1e6;

// EoS ------

// Ideal
double gammaFac = 5.0/3.0;
double mass = 1.67e-27;
double constResis = 3;
double mass_frac_n = 0.99;
double mass_frac_i = 1.0 - mass_frac_n;

// Source update and DC ------

bool doDC = true;
bool doSourceUpdate = true;
int sourceTimeRatio = 1;
int impExRatio = 1;
// vector<implicitSource> implicitSources = {ohmic_diffusion};
vector<implicitSource> implicitSources = {};
vector<SourceFuncEx> exSourceFuncs = {w_evolution_func, heating};

// evolver & BCs ------

BCFunc BC = BohmBCs2;
SLICEvolver evolver(BC, mesh);
shared_ptr<Evolver> evolverPtr = make_shared<SLICEvolver>(evolver);

// recorder and exporter ------

double recordingDelayTime = (tMax-tMin)/20;
Exporter exporter = ExportPIP;

// folders and mixture ------

string mixName = "HFusion";
string resultsFolder = "./output/detachment/";
string dataFolder = "./EoSData/HFusion4_S1/";

// forcing initial time steps ------

/* the step number up until which time step is forced */
int forced_step_lim = 1e6;

/* the ratio to lower time step by */
double forced_ratio = 0.0001;

int main(void){

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
    Vec2D uInit = initDetachment(initParams, mixName, mesh, sysPtr, BC);

    Recorder recorder;
    recorder.setDelayTime(recordingDelayTime);
    shared_ptr<Recorder> recorderPtr = make_shared<Recorder>(recorder);

    Simulation sim(uInit, evolverPtr, sysPtr, recorderPtr, BC, exporter, mesh, resultsFolder, tMax);
    sim.setDoSourceUpdate(doSourceUpdate);
    sim.setSourceTimeRatio(sourceTimeRatio, impExRatio);
    sim.setImplicitSources(implicitSources);
    sim.setExplicitSourceFuncs(exSourceFuncs);
    sim.setDoDC(doDC);

    sim.forceRecordAll();
    sim.enableProgressUpdate(0.1);
    sim.inform();

    cout << "Starting simulation. \n\n" << endl;
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
        cout << "t: " << t << endl;
    }

    cout << "Simulation Completed." << endl;
    sim.forceRecordAll();
    sim.exportAll(resultsFolder);


    return 0;

}

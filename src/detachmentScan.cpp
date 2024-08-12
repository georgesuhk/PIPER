#include "simulation.hpp"
#include "detachInit.hpp"

using namespace std;
const int cellVarsNums = 12;

// grid ------

int nCellsX = 30, nCellsY = 2;
double xMin = 0, xMax = 1;
double yMin = 0, yMax = 0.2;
double tMin = 0, tMax = 1;

double dvar_tol = 0.00025;


Mesh2D mesh(xMin, xMax, nCellsX, yMin, yMax, nCellsY);


// init params ------

const int omp_threads = 1;
double rho_SF = 1.8e-7;
double p_SF = 1.37e-3;
// vector<double> upstream_n_array = {1e19, 1.25e19, 1.5e19, 1.75e19, 2e19, 2.25e19, 2.5e19, 2.75e19, 3e19, 3.25e19, 3.5e19, 3.75e19, 4e19, 4.25e19, 4.5e19, 4.75e19, 5e19, 5.25e19, 5.5e19, 5.75e19, 6e19, 6.25e19, 6.5e19, 6.75e19, 7e19};
// vector<double> upstream_n_array = {1e19, 1.5e19, 2e19, 2.5e19,  3e19, 3.5e19, 4e19, 4.5e19,  5e19, 5.5e19, 6e19,  6.5e19,  7e19, 7.5e19, 8e19, 8.5e19, 9e19, 9.5e19, 10e19};
// vector<double> upstream_n_array = {1e19, 2e19, 3e19, 4e19, 5e19, 6e19,  7e19, 8e19, 9e19, 10e19, 11e19};


vector<double> upstream_n_array = {7e19};

int maxSteps = 30000;

// EoS ------
double mass = 1.67e-27;

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

double recordingDelayStep = 5;
Exporter exporter = ExportPIP;

// folders and mixture ------

string mixName = "HFusion";
string resultsFolder = "./output/detachment/";
string dataFolder = "./EoSData/HFusion4_S4/";

// forcing initial time steps ------

/* the step number up until which time step is forced */
int forced_step_lim = 1e6;

/* the ratio to lower time step by */
double forced_ratio = 0.0001;

int main(void){
    Scalar1D target_T_array, target_flux_array, output_upstream_n_array;

    /* systems that stay constant between sims */

    TabEoS EoSTab;
    EoSTab.genFromData(mesh, {"pressure","densities","T","Cs","e","gamma","resis","thermCon","mass_frac_e","mass_frac_n","mass_frac_i"}, dataFolder, ',');
    shared_ptr<EoS> EoSPtr = make_shared<TabEoS>(EoSTab);

    PIP0_Calcs sysCalcs(EoSPtr);
    shared_ptr<SysCalcs> sysPtr = make_shared<PIP0_Calcs>(sysCalcs);

    Recorder recorder;
    recorder.setDelayStep(recordingDelayStep);
    shared_ptr<Recorder> recorderPtr = make_shared<Recorder>(recorder);

    /* sim dependent systems */
    for (int simNum = 0; simNum < upstream_n_array.size(); simNum++){
        double upstream_n = upstream_n_array[simNum];
        double upstream_rho = upstream_n * mass;
        double upstream_T = sysPtr->interp_T(upstream_rho, p_SF);
        vector<double> initParams = {1*p_SF, upstream_T, 5000, 1*sqrt(p_SF), 0};

        cout << "\n" << endl;
        cout << "upstream n: " << upstream_n << endl;
        cout << "upstream rho: " << upstream_rho << endl;
        cout << "upstream T: " << upstream_T << endl;
        cout << "dvar tol: " << dvar_tol << endl;

        Vec2D uInit = initDetachment(initParams, mixName, mesh, sysPtr, BC);

        Simulation sim(uInit, evolverPtr, sysPtr, recorderPtr, BC, exporter, mesh, resultsFolder, tMax);
        sim.setDoSourceUpdate(doSourceUpdate);
        sim.setSourceTimeRatio(sourceTimeRatio, impExRatio);
        sim.setImplicitSources(implicitSources);
        sim.setExplicitSourceFuncs(exSourceFuncs);
        sim.setDoDC(doDC);

        sim.forceRecordAll();

        cout << "Starting simulation. \n\n" << endl;
        double t = tMin, step = 1;

        double dvar_max = 1e100;
        Vec2D uStart, uEnd;
        int j = 1;
        double dRho, dT, dv = 0;

        while (dvar_max > dvar_tol && step <= maxSteps){
            uStart = sim.get_u();
            
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
            uEnd = sim.get_u();

            dRho = fabs((uEnd[mesh.nCellsX][j][0] - uStart[mesh.nCellsX][j][0]) / uStart[mesh.nCellsX][j][0]);
            dv = fabs((uEnd[mesh.nCellsX][j][1]/uEnd[mesh.nCellsX][j][0] - uStart[mesh.nCellsX][j][1]/uStart[mesh.nCellsX][j][0]) / (uStart[mesh.nCellsX][j][1]/uStart[mesh.nCellsX][j][0]));
            dT = fabs( (sysCalcs.interp_T(uEnd[mesh.nCellsX][j]) - sysCalcs.interp_T(uStart[mesh.nCellsX][j])) / sysCalcs.interp_T(uStart[mesh.nCellsX][j]));
            dvar_max = max({dRho, dv, dT});

            cout << "dRho: " << dRho << endl;
            cout << "dv: " << dv << endl;
            cout << "dT: " << dT << endl;

            cout << "dvar max: " << dvar_max << endl;
        }

        cout << "Simulation Completed." << endl;
        sim.forceRecordAll();
        sim.exportAll(resultsFolder);

        // recording target variables
        CellVec target_u = sim.get_u()[mesh.nCellsX][j];
        target_T_array.push_back(sysCalcs.interp_T(target_u)*kB/elementaryCharge);
        target_flux_array.push_back(sysCalcs.get_vix(target_u) * (target_u[0] / mass));
        output_upstream_n_array.push_back(upstream_n);

        cout << "target_T_array: " << target_T_array << endl;
        cout << "target_flux_array: " << target_flux_array << endl;

        singleColExport(output_upstream_n_array, "upstream_n", resultsFolder);
        singleColExport(target_T_array, "target_T", resultsFolder);
        singleColExport(target_flux_array, "target_flux", resultsFolder);


    }

    return 0;

}

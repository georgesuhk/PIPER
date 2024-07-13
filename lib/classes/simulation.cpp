#include "simulation.hpp"

// ACCESS ======

/* u (TopVec) */
void Simulation::set_u(Vec2D& input_u){
    u = input_u;
}

Vec2D& Simulation::get_u(){
    return u;
}

/* evolver */
void Simulation::setEvolverPtr(shared_ptr<Evolver> inputEvolverPtr){
    evolverPtr = inputEvolverPtr;
}

shared_ptr<Evolver> Simulation::getEvolverPtr(){
    return evolverPtr;
}

/* sysCalc */
void Simulation::setSysPtr(shared_ptr<SysCalcs> inputSysPtr){
    sysPtr = inputSysPtr;
}

shared_ptr<SysCalcs> Simulation::getSysPtr(){
    return sysPtr;
}

/* recorder */
void Simulation::setRecorderPtr(shared_ptr<Recorder> inputRecorderPtr){
    recorderPtr = inputRecorderPtr;
}

shared_ptr<Recorder> Simulation::getRecorderPtr(){
    return recorderPtr;
}

/* BC */
void Simulation::setBC(BCFunc inputBCs){
    BC = inputBCs;
}

/* mesh */
void Simulation::setMesh(Mesh2D inputMesh){
    mesh = inputMesh;
}

Mesh2D Simulation::getMesh(){
    return mesh;
}

/* exporter */
void Simulation::setExporter(Exporter inputExporter){
    exporter = inputExporter;
}

/* simulation time and step */
double Simulation::getTime(){
    return t;
}; 

double Simulation::getStep(){
    return step;
}

void Simulation::setResultsFolder(string inputResultsFolder){
    resultsFolder = inputResultsFolder;
}

string Simulation::getResultsFolder(){
    return resultsFolder;
}



// FUNCTIONALITY ======

double Simulation::get_min_dt(){
    evolverPtr->updateCh(u, sysPtr, mesh);
    return evolverPtr->getTimeStep(mesh);
}

void Simulation::update_dt(){
    dt = get_min_dt();
}

void Simulation::force_set_dt(double input_dt){
    dt = input_dt;
}



void Simulation::evolve(){

    // applying boundary conditions
    BC(u, mesh, sysPtr);


    // getting time step and updating time
    t += dt;
    step += 1;


    // give progress update if giveProgressUpdate enabled
    if (giveProgressUpdate){
        if (t < tMax){
            progressUpdateCounter += dt;
            if (progressUpdateCounter > progressUpdateTime){
                cout << "\rProgress: " << t/tMax * 100 << " %      " << std::flush;
                progressUpdateCounter = 0;
            }
        } else {
            cout << "\rProgress: " << 100 << " %      " << std::endl;
        }
    }
    /* Caching */
    sysPtr->getEoSPtr()->cacheAll(u, mesh);

    // EVOLVE MATERIALS - Strang Splitting ======

    /* source term evolution - S^(t/2) */
    if (doSourceUpdate){
        sourceUpdates(dt/2);
    }

    /* conservative update - C^(t)*/
    evolverPtr->evolveMat(u, sysPtr, mesh, dt, 'x');
    sysPtr->getEoSPtr()->cacheAll(u, mesh);

    evolverPtr->evolveMat(u, sysPtr, mesh, dt, 'y');
    sysPtr->getEoSPtr()->cacheAll(u, mesh);

    /* source term evolution - S^(t/2) */
    if (doSourceUpdate){
        sourceUpdates(dt/2);
    }

    // RECORDING MATERIALS ======
    recorderPtr->update(dt, t, step, u);

    // cout << u[round(mesh.nCellsX/2)][2] << endl;
}

void Simulation::sourceUpdates(double input_dt){

    // implicit source update time step
    double source_dt_imp = input_dt / source_time_ratio;
    double source_dt_ex = source_dt_imp / imp_ex_time_ratio;

    // how many times the implicit update has been applied
    int impSourceUpdateCounter = 0;
    int exSourceUpdateCounter = 0;
    

    // diffusion BC
    string diffBC = "Neumann";

    while (impSourceUpdateCounter < source_time_ratio){
        impSourceUpdateCounter += 1;

        /* implicit updates */
        for (implicitSource& source : implicitSources){
            source(u, sysPtr, mesh, source_dt_imp, diffBC, BC);
        }

        /* explicit updates - can be multiple for each implicit update */
        while (exSourceUpdateCounter < imp_ex_time_ratio){
            exSourceUpdateCounter += 1;
            for (SourceFuncEx& sourceFuncEx : explicitSourceFuncs){
                explicitSolver(u, sysPtr, mesh, source_dt_ex, sourceFuncEx, BC);
            }
        }

    }
}

void Simulation::forceRecordAll(){
    recorderPtr->record(t, step, u);
}

void Simulation::exportAll(string resultsFolder){
    recorderPtr->exportData(exporter, mesh, sysPtr, resultsFolder);
}



// INFORMATION AND DEBUGGING ======

void Simulation::inform(){
    int colwidth = 30;
    cout << "\n-------------" << endl;

    //descriptors for the system
    string sys = sysPtr->getSysName();

    //outputting descrip
    cout << "This is a " << sys << " simulation." << endl;

    if (doDC){
        cout << "\n Divergence Cleaning is enabled. \n " << endl;
    }

    cout << "- " << endl;
    cout << left << setw(colwidth) << "EoS Type: " << sysPtr->getEoSPtr()->getEoSType() << endl;
    cout << left << setw(colwidth) << "EvolverType: " << evolverPtr->getEvolverType() << endl;
    cout << left << setw(colwidth) << "Recording time delay: " << recorderPtr->getDelayTime() << "s" << endl;
    cout << "------------- " << endl;

    cout << "Mesh and time limits:" << endl;

    cout << left << setw(colwidth) << "nCellsX: " << mesh.nCellsX << ", nCellsY: " << mesh.nCellsY << endl;
    cout << left << setw(colwidth) << "dx: " << mesh.dx << ", dy: " << mesh.dy << endl;     
    cout << left << setw(colwidth) << "tMax: " << tMax << endl;

    cout << "------------- \n" << endl;
}

void Simulation::enableProgressUpdate(double updatePercentage){
    giveProgressUpdate = true;
    progressUpdateTime = tMax * updatePercentage;
}

void Simulation::informFinished(){
    if (giveProgressUpdate){
        cout << "\rProgress: 100 %               " << std::flush;
    }

    cout << "\n\nSimulation complete." << endl;
    //could add exporting here
}




// CONTROL ======

void forceTimeStep(double input_dt){

}




// OTHER SYSTEMS ======

// SOURCE UPDATE ------

void Simulation::setDoSourceUpdate(bool inputDoSource){
    doSourceUpdate = inputDoSource;

    if (doSourceUpdate){
        cout << "\n" << endl;
        cout << "Source step updates enabled. With time step ratio: " << source_time_ratio << endl;
        cout << "\n" << endl;
    }
}

bool Simulation::getDoSourceUpdate(){
    return doSourceUpdate;
}

void Simulation::setSourceTimeRatio(int inputSourceRatio){
    source_time_ratio = inputSourceRatio;
}

void Simulation::setSourceTimeRatio(int impSourceRatio, int exSourceRatio){
    source_time_ratio = impSourceRatio;
    imp_ex_time_ratio = exSourceRatio;
}

int Simulation::getSourceTimeRatio(){
    return source_time_ratio;
}

void Simulation::setImpExTimeRatio(int inputTimeRatio){
    imp_ex_time_ratio = inputTimeRatio;
}

int Simulation::getImpExTimeRatio(){
    return source_time_ratio;
}

void Simulation::setImplicitSources(vector<implicitSource> inputSources){
    implicitSources = inputSources;
}

void Simulation::setExplicitSourceFuncs(vector<SourceFuncEx> inputSources){
    explicitSourceFuncs = inputSources;
}

void Simulation::setExplicitSolver(ExplicitSolver inputSolver){
    explicitSolver = inputSolver;
}


vector<implicitSource>& Simulation::getImplicitSources(){
    return implicitSources;
}


// DIVERGENCE CLEANING ------

void Simulation::setDoDC(bool inputDoDC){
    doDC = inputDoDC;
    evolverPtr->setDoDC(inputDoDC);
}

bool Simulation::getDoDC(){
    return doDC;
}


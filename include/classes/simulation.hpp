#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "evolver.hpp"
#include "recorders.hpp"
#include "init.hpp"
#include "sourceUpdates.hpp"

class Simulation {
    public:
        //default constructor
        Simulation(){};

        //full constructor
        Simulation(Vec2D input_u, shared_ptr<Evolver> inputEvolverPtr, shared_ptr<SysCalcs> inputSysPtr, shared_ptr<Recorder> inputRecorder, 
        BCFunc inputBC, Exporter inputExporters, Mesh2D inputMesh, string inputResultsFolder, double input_tMax): 
        u(input_u),
        evolverPtr(inputEvolverPtr),
        sysPtr(inputSysPtr),
        recorderPtr(inputRecorder),
        BC(inputBC),
        exporter(inputExporters),
        mesh(inputMesh),
        resultsFolder(inputResultsFolder),
        tMax(input_tMax){};



        // ACCESS ======

        /* u (TopVec) */
        void set_u(Vec2D& input_u);
        Vec2D& get_u();

        /* evolver */
        void setEvolverPtr(shared_ptr<Evolver> inputEvolverPtr);
        shared_ptr<Evolver> getEvolverPtr();

        /* sysCalc */
        void setSysPtr(shared_ptr<SysCalcs> inputSysPtr);
        shared_ptr<SysCalcs> getSysPtr();
        
        /* recorder */
        void setRecorderPtr(shared_ptr<Recorder> inputRecorderPtr);
        shared_ptr<Recorder> getRecorderPtr();

        /* BC */
        void setBC(BCFunc inputBCs);

        /* mesh */
        void setMesh(Mesh2D inputMesh);
        Mesh2D getMesh();

        /* exporter */
        void setExporter(Exporter inputExporter);

        /* simulation time and steps */
        double getTime();
        double getStep();

        /* results folder */
        void setResultsFolder(string inputResultsFolder);
        string getResultsFolder();



        // FUNCTIONALITY ======

        /* calculates and returns the minimum dt across all the evolvers in the simulation */
        double get_min_dt();

        /* simply returns the dt value stored in the sim */
        double get_dt();

        /* updates the stored dt with min dt */
        void update_dt();

        /* forces dt to a particular value */
        void force_set_dt(double input_dt);

        /* applies conservative evolution to all cells in the simulation */
        void evolve();

        /* applies source terms evolution to all cells in the simulation */
        void sourceUpdates(double dt);

        /**
         * forces all the recorders to take a snapshot of the system at the current time step
         * By passes the delay within the recorders
         */
        void forceRecordAll();

        /* export all recorded states */
        void exportAll(string resultsFolder);

        void setDoInSimExport(bool input);
        void setExportGap(double input);



        // INFORMATION AND DEBUGGING ======

        /* gives a brief of the simulation being ran */
        void inform();

        /* enables the printing of the simulation progress */
        void enableProgressUpdate(double updatePercentage);
        
        void informFinished();



        // OTHER SYSTEMS ======

        // SOURCE UPDATE ------

        /* toggling source update */
        void setDoSourceUpdate(bool inputDoSource);
        bool getDoSourceUpdate();

        /* setting the ratio been source update and conservative update */
        /* default is 1:1 */
        void setSourceTimeRatio(int inputSourceRatio);

        /* first input is the ratio between source update and conserv update, second input the ratio between
        implicit and explicit update */
        void setSourceTimeRatio(int impSourceRatio, int exSourceRatio);

        int getSourceTimeRatio();

        /* setting the ratio been implicit update and explicit update */
        /* default is 1:1 */
        void setImpExTimeRatio(int inputTimeRatio);
        int getImpExTimeRatio();


        // ADDING SOURCES ------

        /* adding implicit source updates */
        void setImplicitSources(vector<implicitSource> inputSources);
        vector<implicitSource>& getImplicitSources();

        /* adding explicit source update functions */
        void setExplicitSourceFuncs(vector<SourceFuncEx> inputSources);

        /* sets which solver algorithm to use for the explicit solver */
        void setExplicitSolver(ExplicitSolver inputSolver);


        // DIVERGENCE CLEANING ------

        /* toggling div clean */
        void setDoDC(bool inputDoDC);
        bool getDoDC();

    protected:

        // Core components ------

        Vec2D u;
        shared_ptr<Evolver> evolverPtr;
        shared_ptr<SysCalcs> sysPtr;
        shared_ptr<Recorder> recorderPtr;
        BCFunc BC;
        Exporter exporter;
        Mesh2D mesh;


        // source term update ------

        bool doSourceUpdate = false;
        int source_time_ratio = 1;
        int imp_ex_time_ratio = 1;

        vector<implicitSource> implicitSources;
        vector<SourceFuncEx> explicitSourceFuncs;
        ExplicitSolver explicitSolver = RK2;

        

        // control variables ------

        string resultsFolder;
        double dt;
        double t = 0;
        int step = 0;
        double tMax;
        bool doDC = false;


        //info and system properties ------
        double giveProgressUpdate = false;
        double progressUpdateCounter = 0;
        double progressUpdateTime;

        bool doInSimExport = false;
        double exportCounter = 0;
        double exportTimeGap = 1e50;
};

#endif

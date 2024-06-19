#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "evolver.hpp"
#include "recorders.hpp"
#include "init.hpp"

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

        /* returns the minimum dt across all the evolvers in the simulation */
        double get_min_dt();

        /* evolvers all the TopVec in the simulation */
        void evolve();

        /**
         * forces all the recorders to take a snapshot of the system at the current time step
         * By passes the delay within the recorders
         */
        void forceRecordAll();

        /* export all recorded states */
        void exportAll(string resultsFolder);



        // INFORMATION AND DEBUGGING ======

        /* gives a brief of the simulation being ran */
        void inform();

        /* enables the printing of the simulation progress */
        void enableProgressUpdate(double updatePercentage);
        
        void informFinished();


        // CONTROL ======

        /* toggling source update */

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
        

        // control variables ------

        string resultsFolder;
        double t = 0;
        int step = 0;
        double tMax;
        bool doDC = false;


        //info and system properties ------
        double giveProgressUpdate = false;
        double progressUpdateCounter = 0;
        double progressUpdateTime;
};

#endif

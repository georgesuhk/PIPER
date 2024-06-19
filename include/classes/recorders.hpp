#ifndef RECORDERS_HPP
#define RECORDERS_HPP

#include "evolver.hpp"
#include "exporters.hpp"

class Recorder {
    public:
        Recorder(): delayTime(numeric_limits<double>::quiet_NaN()), delayStep(-1){};

        /* sets the time interval with which recordings are taken */
        void setDelayTime(double delay);

        /* set the step interval with which recordings are taken */
        void setDelayStep(int delay);

        /* returns recording time interval */
        double getDelayTime();

        /* returns recording step interval */
        double getDelayStep();

        /* return the steps for which recording happened */
        vector<int> getRecordSteps();

        /* returns the times at which recording happened */
        vector<double> getRecordTimes();

        /* returns the stored time series of TopVec */
        vector<Vec2D> getStoredVec();

        //updating internal time/step, to check if a snapshot needs to be taken
        void update(double dt, double t, int step, Vec2D& u);

        //recording data
        void record(double t, int step, Vec2D& u);

        //exporting data
        void exportData(Exporter exporter, Mesh2D mesh, shared_ptr<SysCalcs> sysPtr, string folder);
        

    protected:
        // the counters which manage if a new snapshot needs to be taken
        double delayCounterTime;
        int delayCounterStep;

        // the delay in time/steps between recordings
        double delayTime;
        int delayStep;

        // list of all the times/steps at which recordings were taken
        vector<double> recordTimes; 
        vector<int> recordSteps;

        // stored TopVec of variables
        vector<Vec2D> storedVec; //stores all top vecs for 2D sims

};

#endif
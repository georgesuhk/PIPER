#include "recorders.hpp"

void Recorder::setDelayTime(double delay){
    delayTime = delay;
    delayCounterTime = delayTime;
}

void Recorder::setDelayStep(int delay){
    delayStep = delay;
    delayCounterStep = delayStep;
}

double Recorder::getDelayTime(){
    return delayTime;
}

double Recorder::getDelayStep(){
    return delayStep;
}

vector<int> Recorder::getRecordSteps(){
    return recordSteps;
}

vector<double> Recorder::getRecordTimes(){
    return recordTimes;
}

vector<Vec2D> Recorder::getStoredVec(){
    return storedVec;
}

void Recorder::update(double dt, double t, int step, Vec2D& u){
    
    if (!isnan(delayTime)){
        delayCounterTime -= dt; //advance time for delay counter

        if (delayCounterTime <= 0){
            record(t, step, u);
            delayCounterTime = delayTime; //reset counter values
            delayCounterStep = delayStep;
            return;  //early termination so steps counter no need to trigger
            
        }
    } 
    if (delayStep != -1){
        delayCounterStep -= 1;         //advance step for delay counter
        if (delayCounterStep <= 0){
            record(t, step, u);
            delayCounterTime = delayTime; //reset counter values
            delayCounterStep = delayStep;
            return;
        }
    }
}

void Recorder::record(double t, int step, Vec2D& u){
    storedVec.push_back(u);
    recordTimes.push_back(t);
    recordSteps.push_back(step);
}

void Recorder::exportData(Exporter exporter, Mesh2D mesh, shared_ptr<SysCalcs> sysPtr, string folder){
    exporter(storedVec, recordTimes, mesh, sysPtr, folder);
}

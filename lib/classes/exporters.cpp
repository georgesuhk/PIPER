#include "exporters.hpp"

// FUNCTIONS USED IN EXPORTING ======

double getTotalU(Vec2D& u, Mesh2D& mesh){
    double U = 0;
    for (int i = 1; i < mesh.nCellsX+1; i++){
        for (int j = 1; j < mesh.nCellsY+1; j++){
            U += u[i][j][4];
        }
    }    
    return U;
}



// EXPORTERS ======
void ExportFIP(vector<Vec2D>& uTimeSeries, vector<double>& recordedTimes, Mesh2D mesh, shared_ptr<SysCalcs> sysPtr, string folder){

    cout << "Exporting using 'ExportDefault' to folder: " << folder << endl;

    // cell specific variables
    ofstream rhoFile(folder+"rho.csv");
    ofstream vxFile(folder+"vx.csv");
    ofstream vyFile(folder+"vy.csv");
    ofstream vzFile(folder+"vz.csv");
    ofstream pFile(folder+"p.csv");
    ofstream eFile(folder+"e.csv");
    ofstream TFile(folder+"T.csv"); 
    ofstream BxFile(folder+"Bx.csv");
    ofstream ByFile(folder+"By.csv");
    ofstream BzFile(folder+"Bz.csv");

    // total variables
    ofstream totalUFile(folder+"totalU.csv");
    
    // grid
    ofstream gridXFile(folder+"gridX.csv");
    ofstream gridYFile(folder+"gridY.csv");

    // recording times
    ofstream recordTimeFile(folder+"recordTimes.csv");

    // WRITING ======

    // non-cell specific variables ------
    
    totalUFile << "[";
    recordTimeFile << "[";

    for (int timeIdx = 0; timeIdx < int(uTimeSeries.size()); timeIdx++){
        if (timeIdx != int(uTimeSeries.size())-1){
            totalUFile << getTotalU(uTimeSeries[timeIdx], mesh) << ",";
            recordTimeFile << recordedTimes[timeIdx] << ",";
        } else {
            totalUFile << getTotalU(uTimeSeries[timeIdx], mesh) << "]" << endl;
            recordTimeFile << recordedTimes[timeIdx] << "]" << endl;
        }
    } 


    // cell specific variables ------
    double x, y;

    //recording x and y grids
    for (int j = 0; j < mesh.nCellsY; j++){
        for (int i = 0; i < mesh.nCellsX; i++){
            x = mesh.xMin + i * mesh.dx;
            y = mesh.yMin + j * mesh.dy;

            if (i != mesh.nCellsX-1){
                gridXFile << x << ",";
                gridYFile << y << ",";
            } else { //final value
                gridXFile << x << endl;
                gridYFile << y << endl;
            }
        }
    }

    double rho, vx, vy, vz, p, e, T, Bx, By, Bz; 
    for (int j = 1; j < mesh.nCellsY+1; j++){
        for (int i = 1; i < mesh.nCellsX+1; i++){
        
            rhoFile << "[";
            vxFile << "[";
            vyFile << "[";
            vzFile << "[";
            pFile << "[";
            eFile << "[";
            TFile << "[";
            BxFile << "[";
            ByFile << "[";
            BzFile << "[";

            Vec2D u;
            CellVec uPrim;
            for (int timeIdx = 0; timeIdx < int(uTimeSeries.size()); timeIdx++){
                u = uTimeSeries[timeIdx];

                uPrim = sysPtr->conservToPrim(u[i][j]);

                rho = uPrim[0];
                vx = uPrim[1];
                vy = uPrim[2];
                vz = uPrim[3];
                p = uPrim[4];
                e = sysPtr->get_e_Prim(uPrim);
                T = sysPtr->get_T_Prim(uPrim);
                Bx = uPrim[5];
                By = uPrim[6];
                Bz = uPrim[7];

                if (timeIdx != int(uTimeSeries.size())-1){
                    rhoFile << rho << ",";
                    vxFile << vx << ",";
                    vyFile << vy << ",";
                    vzFile << vz << ",";
                    pFile << p << ",";
                    eFile << e << ",";
                    TFile << T << ",";
                    BxFile << Bx << ",";
                    ByFile << By << ",";
                    BzFile << Bz << ",";

                }
                else {
                    if (i != mesh.nCellsX){
                        rhoFile << rho << "],";
                        vxFile << vx << "],";
                        vyFile << vy << "],";
                        vzFile << vz << "],";
                        pFile << p << "],";
                        eFile << e << "],";
                        TFile << T << "],";
                        BxFile << Bx << "],";
                        ByFile << By << "],";
                        BzFile << Bz << "],"; 
                    }
                    else {
                        rhoFile << rho << "]" << endl;
                        vxFile << vx << "]" << endl;
                        vyFile << vy << "]" << endl;
                        vzFile << vz << "]" << endl;
                        pFile << p << "]" << endl;
                        eFile << e << "]" << endl;
                        TFile << T << "]" << endl;
                        BxFile << Bx << "]" << endl;
                        ByFile << By << "]" << endl;
                        BzFile << Bz << "]" << endl;
                    }
                }
            }
        }
    }

    cout << "Exporting finished." << endl;

}

void ExportPIP(vector<Vec2D>& uTimeSeries, vector<double>& recordedTimes, Mesh2D mesh, shared_ptr<SysCalcs> sysPtr, string folder){

    cout << "Exporting using 'ExportDefault' to folder: " << folder << endl;

    // cell specific variables
    ofstream rhoFile(folder+"rho.csv");
    ofstream vxFile(folder+"vx.csv");
    ofstream vyFile(folder+"vy.csv");
    ofstream vzFile(folder+"vz.csv");
    ofstream pFile(folder+"p.csv");
    ofstream eFile(folder+"e.csv");
    ofstream TFile(folder+"T.csv"); 
    ofstream BxFile(folder+"Bx.csv");
    ofstream ByFile(folder+"By.csv");
    ofstream BzFile(folder+"Bz.csv");
    ofstream wxFile(folder+"wx.csv");
    ofstream wyFile(folder+"wy.csv");
    ofstream wzFile(folder+"wz.csv");
    ofstream mass_frac_n_File(folder+"mass_frac_n.csv");
    ofstream mass_frac_i_File(folder+"mass_frac_i.csv");
    ofstream vnxFile(folder+"vnx.csv");
    ofstream vixFile(folder+"vix.csv");
    ofstream gammaFile(folder+"gamma.csv");
    ofstream resisFile(folder+"resis.csv");


    // total variables
    ofstream totalUFile(folder+"totalU.csv");
    
    // grid
    ofstream gridXFile(folder+"gridX.csv");
    ofstream gridYFile(folder+"gridY.csv");

    // recording times
    ofstream recordTimeFile(folder+"recordTimes.csv");

    // WRITING ======

    // non-cell specific variables ------
    
    totalUFile << "[";
    recordTimeFile << "[";

    for (int timeIdx = 0; timeIdx < int(uTimeSeries.size()); timeIdx++){
        if (timeIdx != int(uTimeSeries.size())-1){
            totalUFile << getTotalU(uTimeSeries[timeIdx], mesh) << ",";
            recordTimeFile << recordedTimes[timeIdx] << ",";
        } else {
            totalUFile << getTotalU(uTimeSeries[timeIdx], mesh) << "]" << endl;
            recordTimeFile << recordedTimes[timeIdx] << "]" << endl;
        }
    } 


    // cell specific variables ------
    double x, y;

    //recording x and y grids
    for (int j = 0; j < mesh.nCellsY; j++){
        for (int i = 0; i < mesh.nCellsX; i++){
            x = mesh.xMin + i * mesh.dx;
            y = mesh.yMin + j * mesh.dy;

            if (i != mesh.nCellsX-1){
                gridXFile << x << ",";
                gridYFile << y << ",";
            } else { //final value
                gridXFile << x << endl;
                gridYFile << y << endl;
            }
        }
    }

    double rho, vx, vy, vz, p, e, T, Bx, By, Bz, wx, wy, wz, mass_frac_n, mass_frac_i, vnx, vix, gamma, resis; 
    for (int j = 1; j < mesh.nCellsY+1; j++){
        for (int i = 1; i < mesh.nCellsX+1; i++){
        
            rhoFile << "[";
            vxFile << "[";
            vyFile << "[";
            vzFile << "[";
            pFile << "[";
            eFile << "[";
            TFile << "[";
            BxFile << "[";
            ByFile << "[";
            BzFile << "[";
            wxFile << "[";
            wyFile << "[";
            wzFile << "[";
            mass_frac_n_File << "[";
            mass_frac_i_File << "[";
            vnxFile << "[";
            vixFile << "[";
            gammaFile << "[";
            resisFile << "[";


            Vec2D u;
            CellVec uPrim;
            for (int timeIdx = 0; timeIdx < int(uTimeSeries.size()); timeIdx++){
                u = uTimeSeries[timeIdx];

                uPrim = sysPtr->conservToPrim(u[i][j]);

                rho = uPrim[0];
                vx = uPrim[1];
                vy = uPrim[2];
                vz = uPrim[3];
                p = uPrim[4];
                e = sysPtr->get_e_Prim(uPrim);
                T = sysPtr->get_T_Prim(uPrim);
                Bx = uPrim[5];
                By = uPrim[6];
                Bz = uPrim[7];
                //psi
                wx = uPrim[9];
                wy = uPrim[10];
                wz = uPrim[11];
                mass_frac_i = sysPtr->interp_mass_frac_i_Prim(uPrim);
                mass_frac_n = 1 - mass_frac_i;
                vix = vx + mass_frac_n * wx;
                vnx = vx - mass_frac_i * wx;
                gamma = sysPtr->interp_gamma(u[i][j]);
                resis = sysPtr->interp_Resis(rho, p);

                if (timeIdx != int(uTimeSeries.size())-1){
                    rhoFile << rho << ",";
                    vxFile << vx << ",";
                    vyFile << vy << ",";
                    vzFile << vz << ",";
                    pFile << p << ",";
                    eFile << e << ",";
                    TFile << T << ",";
                    BxFile << Bx << ",";
                    ByFile << By << ",";
                    BzFile << Bz << ",";
                    wxFile << wx << ",";
                    wyFile << wy << ",";
                    wzFile << wz << ",";
                    mass_frac_n_File << mass_frac_n << ",";
                    mass_frac_i_File << mass_frac_i << ",";
                    vnxFile << vnx << ",";
                    vixFile << vix << ",";
                    gammaFile << gamma << ",";
                    resisFile << resis << ",";


                }
                else {
                    if (i != mesh.nCellsX){
                        rhoFile << rho << "],";
                        vxFile << vx << "],";
                        vyFile << vy << "],";
                        vzFile << vz << "],";
                        pFile << p << "],";
                        eFile << e << "],";
                        TFile << T << "],";
                        BxFile << Bx << "],";
                        ByFile << By << "],";
                        BzFile << Bz << "],"; 
                        wxFile << wx << "],";
                        wyFile << wy << "],";
                        wzFile << wz << "],";
                        mass_frac_n_File << mass_frac_n << "],";
                        mass_frac_i_File << mass_frac_i << "],";
                        vnxFile << vnx << "],";
                        vixFile << vix << "],";
                        gammaFile << gamma << "],";
                        resisFile << resis << "],";

                    }
                    else {
                        rhoFile << rho << "]" << endl;
                        vxFile << vx << "]" << endl;
                        vyFile << vy << "]" << endl;
                        vzFile << vz << "]" << endl;
                        pFile << p << "]" << endl;
                        eFile << e << "]" << endl;
                        TFile << T << "]" << endl;
                        BxFile << Bx << "]" << endl;
                        ByFile << By << "]" << endl;
                        BzFile << Bz << "]" << endl;
                        wxFile << wx << "]" << endl;
                        wyFile << wy << "]" << endl;
                        wzFile << wz << "]" << endl;
                        mass_frac_n_File << mass_frac_n << "]" << endl;
                        mass_frac_i_File << mass_frac_i << "]" << endl;
                        vnxFile << vnx << "]" << endl;
                        vixFile << vix << "]" << endl;
                        gammaFile << gamma << "]" << endl;
                        resisFile << resis << "]" << endl;

                    }
                }
            }
        }
    }

    cout << "Exporting finished." << endl;

}


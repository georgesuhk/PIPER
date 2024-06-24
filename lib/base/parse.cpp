#include "parse.hpp"

// PARSE ======

vector<double> singleColParse(string filename){
    vector<double> data;
    ifstream file(filename);

    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
    }

    string line;
    while (getline(file, line)){
        data.push_back(stod(line));
    }

    return data;
}



// EXPORT ======

void singleColExport(vector<double> data, string dataName, string folder){
    ofstream dataFile(folder+dataName+".csv");

    for (int i = 0; i < int(data.size()); i++){
        if (i != int(data.size())-1){
            dataFile << data[i] << " ";
        } else {
            dataFile << data[i] << endl;
        }
    } 
}

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

void parseLine(const std::string& line, std::vector<std::string>& parsedData) {
    std::stringstream ss(line);
    std::string item;
    
    while (std::getline(ss, item, ',')) {
        parsedData.push_back(item);
    }
}


Scalar1D easyParseArray(string filename, char delimiter){
    Scalar1D array;
    ifstream file(filename);

    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
    }

    string line;
    getline(file, line);
    stringstream ss(line);
    string item;
    
    while (getline(ss, item, ',')) {
        array.push_back(stod(item));
    }
     
    return array;
}

Scalar2D easyParseTable(string filename, char delimiter){
    vector<vector<double>> tabularData;
    ifstream file(filename);

    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
    }

    string line;
    
    while (getline(file, line)) {
        vector<string> row;
        stringstream ss(line);
        string cell;

        vector<double> innerVec;
        while (getline(ss, cell, delimiter)) {
            innerVec.push_back(stod(cell));
        }
        tabularData.push_back(innerVec);        
    }

    return tabularData;
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

void easyExport(Scalar1D& data, string dataName, string folder){
    ofstream dataFile(folder+dataName+".csv");

    for (int i = 0; i < data.size(); i++){
        if (i != data.size()-1){
            dataFile << data[i] << ",";
        } else { //final value
            dataFile << data[i] << endl;
        }
    }
}

void easyExport(Scalar2D& data, string dataName, string folder){
    ofstream dataFile(folder+dataName+".csv");

    for (int j = 0; j < data.size(); j++){
        for (int i = 0; i < data[0].size(); i++){
            if (i != data[0].size()-1){
                dataFile << data[j][i] << ",";
            } else { //final value
                dataFile << data[j][i] << endl;
            }
        }
    }
}

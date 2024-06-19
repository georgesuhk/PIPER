#include "opOverload.hpp"


//printing to cout for CellVec
ostream& operator<<(ostream& os, CellVec cellVec){
    os << "[";

    //for each variable
    for (double& value : cellVec){
        os << value << " ";
    }
    os << "]" << endl;

    return os;
}

//printing to cout for vec1D
ostream& operator<<(ostream& os, Vec1D vec1D){

    for(CellVec& vector : vec1D){

        os << "[";

        //for each variable
        for (double& value : vector){
            os << value << " ";
        }
        os << "]" << endl;
    }

    os << endl;
    return os;
}


//printing to cout for TopVec2D
ostream& operator<<(ostream& os, Vec2D vec2D){

    //for each x value
    for (vector<CellVec>& vec1D : vec2D){

        //for each y value 
        for(CellVec& cellVec : vec1D){
            os << "[";

            //for each variable
            for (double& value : cellVec){
                os << value << " ";
            }
            os << "]" << endl;
        }
        os << endl;
    }
    return os;
}

//printing for intArray
ostream& operator<<(ostream& os, vector<int> intArray){

    os << "[";

    //for each variable
    for (int& value : intArray){
        os << value << " ";
    }
    os << "]" << endl;


    return os;
}


ostream& operator<<(ostream& os, array<int, 2> intArray){

    os << "[";

    //for each variable
    for (int& value : intArray){
        os << value << " ";
    }
    os << "]" << endl;

    

    return os;
}


ostream& operator<<(ostream& os, array<double, 2> array){

    os << "[";

    //for each variable
    for (double& value : array){
        os << value << " ";
    }
    os << "]" << endl;

    

    return os;
}
#ifndef OPOVERLOAD_HPP
#define OPOVERLOAD_HPP

#include "settings.hpp"

using namespace std;

// Overload the << operator for a vector of vectors
ostream& operator<<(ostream& os, CellVec cellVec);
ostream& operator<<(ostream& os, Vec1D vec1D);
ostream& operator<<(ostream& os, Vec2D vec2D);


//overloading printing for scalars arrays
ostream& operator<<(ostream& os, vector<int> intArray);
ostream& operator<<(ostream& os, array<int, 2> intArray);
ostream& operator<<(ostream& os, array<double, 2> array);



#endif
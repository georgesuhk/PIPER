#ifndef PARSE_HPP
#define PARSE_HPP

#include "utils.hpp"


/* Contains simple parsing and exporting functions (not to be confused with exporters.hpp) */

// PARSE ======

/* parses data stored in a single column format */
vector<double> singleColParse(string filename);


// EXPORT ======

/* single column exporters for vector<double> data */
void singleColExport(vector<double> data, string dataName, string folder);


#endif 
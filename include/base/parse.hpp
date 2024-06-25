#ifndef PARSE_HPP
#define PARSE_HPP

#include "utils.hpp"


/* Contains simple parsing and exporting functions (not to be confused with exporters.hpp) */

// PARSE ======

/* parses data stored in a single column format */
vector<double> singleColParse(string filename);

/* parser matching to easyExport, for single line data */
Scalar1D easyParseArray(string filename, char delimiter);

/* parser matching to easyExport, for tabular data */
Scalar2D easyParseTable(string filename, char delimiter);

// EXPORT ======

/* single column exporters for vector<double> data */
void singleColExport(vector<double> data, string dataName, string folder);

/* exports as a simple table, delimited by "," */
void easyExport(Scalar1D& data, string dataName, string folder);
void easyExport(Scalar2D& data, string dataName, string folder);

#endif 
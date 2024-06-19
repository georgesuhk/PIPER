#ifndef LIMITER_HPP
#define LIMITER_HPP

#include "mesh.hpp"
#include "utils.hpp"

typedef double (*Limiter)(double);

/* Minbee limiter function */
double minbee(double r);

/* van-Leer limiter function */
double vanLeer(double r);

/* outer function for limiters that returns a top level vector with limiter values for each cell */
Vec2D getLimiterVec(const Vec2D& u, Mesh2D& mesh, Limiter limFunc, char axis);

#endif
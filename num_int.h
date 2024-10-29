#ifndef NUMINT
#define NUMINT

#include "Aux.h"
#include <random>

double simpson1(double (*f)(double) , const double minCoordinate, const double maxCoordinate, const double precision, std::vector<double> &fvals, std::vector<double> &xvals,  long unsigned int &counter);

#endif
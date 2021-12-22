#include <cstdlib>
#include "monte_carlo.h"

// Some functions to draw random numbers which are absent in pure C++

double uniform (double min, double max) {
        double res;
        res = (float) rand()/RAND_MAX;
        res *= (max - min);
        res += min;
        return res;
}

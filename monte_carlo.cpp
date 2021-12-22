#include <iostream>
#include <cstdlib>
#include <ctime>
#include "monte_carlo.h"

using namespace std;

// Some functions to draw random numbers which are absent in pure C++

double uniform (double min, double max) {
        double res;
	res = (double) rand()/ double (RAND_MAX);
        res *= (max - min);
        res += min;
        return res;
}

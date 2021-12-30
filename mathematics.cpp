#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

double interp1d (double x, int n, double * list_x, double * list_y) {

	double a, b, res, t;
	int pos;

	vector<double>::iterator low;

	vector <double> v (list_x, list_x + n);

	low = lower_bound (v.begin(), v.end(), x); // find first element which is above x

	pos = (low - v.begin()) - 1; // return to the previous element


	if (list_x[pos] == x) 
		return list_y[pos]; // check if we already have the value for this element at the grid
	if (pos == n-1) 
		return list_y[n-1]; // check if we are using the last element on the grid

	if (x > list_x[n-1])
		return NAN;         // check if we above the interpolation range
	if (x < list_x[0])
		return NAN;         // check if we below the interpolation range


	t = (x - list_x[pos]) / ( list_x[pos+1] - list_x[pos] );
	a = list_y[pos];
	b = list_y[pos + 1];
	res = a + t * (b - a);
	return res;
}



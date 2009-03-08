#ifndef STIFF_H
#define STIFF_H

#include "matrix.h"
#include <vector>
using namespace std;

extern vector< vector<double> > ode_dat; 
extern int ode_nstep;

#define YMIN 1e-3

void euler_implicit_step(vector<double> y, vector<double> dydx, int dim, double x, double dx, double *yout, 
	double *yerr, void (*derivs)(double, const vector<double> &, vector<double> &, int),
	void (*jacobn)(double, const vector<double> &, double**, int)	);

void rosenbrock_step(vector<double> y, vector<double> dydx, int dim, double x, double dx, vector<double> &yout, 
	vector<double> &yerr, void (*derivs)(double, const vector<double> &, vector<double> &, int),
	void (*jacobn)(double, const vector<double> &, double**, int)	);

void euler_3step(vector<double> y, vector<double> dydx, int dim, double x, double dx, vector<double> &yout, 
	vector<double> &yerr, void (*derivs)(double, const vector<double> &, vector<double> &, int),
	void (*jacobn)(double, const vector<double> &, double**, int)	);

void euler_step(vector<double> y, vector<double> dydx, int dim, double x, double dx, double *yout, 
	double *yerr, void (*derivs)(double, const vector<double> &, vector<double> &, int),
	void (*jacobn)(double, const vector<double> &, double**, int)	);

void stiffint(vector<double> Y0, double x0, double x_max, double eps, double hi,
	void (*derivs)(double, const vector<double> &, vector<double> &, int),
	void (*jacobn)(double, const vector<double> &, double**, int),
	void (*step)(vector<double>, vector<double>, int, double, double, double*, double*,
	    void(*)(double, const vector<double> &, vector<double> &, int),
	    void(*)(double, const vector<double> &, double**, int)),
	int order);

void stiffint_adpt(vector<double> Y0, double x0, double x_max, double eps, double hi, int &nok, int &nbad,
	void (*derivs)(double, const vector<double> &, vector<double> &, int),
	void (*jacobn)(double, const vector<double> &, double**, int),
	void (*step)(vector<double>, vector<double>, int, double, double, vector<double> &, vector<double> &,
	    void(*)(double, const vector<double> &, vector<double> &, int),
	    void(*)(double, const vector<double> &, double**, int)),
	int order);

void odeint(double *Y0, int dim, double x0, double x_max, double eps, double hi, int *nok, int *nbad,
	void (*derivs)(double, double*, double*),
	void (*step)(double*, double*, int, double, double, double*, double*, 
	    void (*)(double, double*, double*)), int order);

void ode_init(int dim, int maxlen);

#endif

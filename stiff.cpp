#include "matrix.h"
#include "stiff.h"
#include <math.h>

int ode_maxlen;

t_numatrix<double> ode_tmp;
vector< vector<double> > ode_dat;
t_LUmatrix<double> ode_jac;
int ode_nstep;

/*
template<class T>
T max(T a, T b)
{
    return A > B ? A : B ;
}

template<class T>
T min(T a, T b)
{
    return A < B ? A : B ;
}
*/


void rosenbrock_step(vector<double> y, vector<double> dydx, int dim, double x, double dx, vector<double> &yout, 
	vector<double> &yerr, void (*derivs)(double, const vector<double> &, vector<double> &, int),
	void (*jacobn)(double, const vector<double> &, double**, int)	)
{
    int i,j;
    vector<double> err(dim);
    vector<double> g1(dim), g2(dim), g3(dim), g4(dim);
    vector<double> ysav = y;
    double xsav = x;
    const double GAM = (1.0/2.0);
    const double A21 = 2.0;
    const double A31 = (48.0/25.0);
    const double A32 = (6.0/25.0);
    const double C21 = -8.0;
    const double C31 = (372.0/25.0);
    const double C32 = (12.0/5.0);
    const double C41 = (-112.0/125.0);
    const double C42 = (-54.0/125.0);
    const double C43 = (-2.0/5.0);
    const double B1 = (19.0/9.0);
    const double B2 = (1.0/2.0);
    const double B3 = (25.0/108.0);
    const double B4 = (125.0/108.0);
    const double E1 = (17.0/54.0);
    const double E2 = (7.0/36.0);
    const double E3 = 0.0;
    const double E4 = (125.0/108.0);
    //const double C1X = (1.0/2.0);
    //const double C2X = (-3.0/2.0);
    //const double C3X = (121.0/50.0);
    //const double C4X = (29.0/250.0);
    const double A2X = 1.0;
    //const double A3X = (3.0/5.0);
    
    jacobn(x, y, ode_jac, dim);
    for(i=0; i<dim; i++)
	for(j=0; j<dim; j++)
	    ode_jac[i][j] *= -1.0;
    for(i=0; i<dim; i++)
	ode_jac[i][i] += 1.0/(GAM*dx);
    ode_jac.decompose();
    

    for(i=0; i<dim; i++)
	g1[i] = dydx[i];
    ode_jac.LUbksb(g1);

    for(i=0; i<dim; i++)
	y[i] = ysav[i] + A21*g1[i];
    
    x = xsav + A2X*dx;
    derivs(x, y, dydx, dim);

    for(i=0; i<dim; i++)
	g2[i] = dydx[i] + C21*g1[i]/dx;
    ode_jac.LUbksb(g2);

    for(i=0; i<dim; i++)
	y[i] = ysav[i] + A31*g1[i] + A32*g2[i];

    x = xsav + A2X*dx;
    derivs(x, y, dydx, dim);

    for(i=0; i<dim; i++)
	g3[i]=dydx[i] + ( C31*g1[i] + C32*g2[i] )/dx;
    ode_jac.LUbksb(g3);
    
    for(i=0; i<dim; i++)
	g4[i] = dydx[i] + ( C41*g1[i] + C42*g2[i] + C43*g3[i])/dx;
    ode_jac.LUbksb(g4);
    
    for(i=0; i<dim; i++)
    {
	ode_tmp[i][0] = yout[i] = ysav[i] + B1*g1[i] + B2*g2[i] + B3*g3[i] + B4*g4[i];
	yerr[i] = E1*g1[i] + E2*g2[i] + E3*g3[i] + E4*g4[i];
    }

}


void euler_implicit_step(vector<double> y, vector<double> dydx, int dim, double x, double dx, double *yout, 
	double *yerr, void (*derivs)(double, const vector<double> &, vector<double> &, int),
	void (*jacobn)(double, const vector<double> &, double**, int)	)
{
    int i,j;
    jacobn(x, y, ode_jac, dim);
    for(i=0; i<dim; i++)
	for(j=0; j<dim; j++)
	    ode_jac[i][j] *= -dx;
    for(i=0; i<dim; i++)
	ode_jac[i][i] += 1.0;

    ode_jac.decompose();
    ode_jac.LUbksb(dydx);
    
    for(i=0; i<dim; i++)
	ode_tmp[i][0] = yout[i] = y[i] + dydx[i]*dx;
}

void euler_step(vector<double> y, vector<double> dydx, int dim, double x, double dx, double *yout, 
	double *yerr, void (*derivs)(double, const vector<double> &, vector<double> &, int),
	void (*jacobn)(double, const vector<double> &, double**, int)	)
{
    int i;
    for(i=0; i<dim; i++)
	ode_tmp[i][0] = yout[i] = y[i] + dydx[i]*dx;
}
void euler_3step(vector<double> y, vector<double> dydx, int dim, double x, double dx, vector<double> &yout, 
	vector<double> &yerr, void (*derivs)(double, const vector<double> &, vector<double> &, int),
	void (*jacobn)(double, const vector<double> &, double**, int)	)
{
    vector<double> C1(dim);
    int i;
    for(i=0; i<dim; i++)
	ode_tmp[i][0] = yout[i] = y[i] + dydx[i]*dx*0.5;
    derivs(x+dx*0.5,yout,C1,dim);
    for(i=0; i<dim; i++)
	ode_tmp[i][1] = yout[i] = yout[i] + C1[i]*dx*0.5;

    for(i=0; i<dim; i++)
    {
	C1[i] = y[i] + dydx[i]*dx;
	yerr[i] = fabs(C1[i]-yout[i]);
    }
}


void stiffint(vector<double> Y0, double x0, double x_max, double eps, double hi,
	void (*derivs)(double,const vector<double> &, vector<double> &, int),
	void (*jacobn)(double,const vector<double> &, double**, int),
	void (*step)(vector<double>, vector<double>, int, double, double, double*, double*,
	    void(*)(double, const vector<double> &, vector<double> &, int),
	    void(*)(double, const vector<double> &, double**, int)),
	int order)
{
    int i;
    double x=x0;
    double h=hi;
    int dim = Y0.size();

    vector<double> Y = Y0;
    vector<double> dY(dim);
    vector<double> tmp_vec(dim);

    ode_nstep = 1;
    
    for(i=0; i<dim; i++)
	ode_dat[i].push_back( Y0[i] );

    ode_dat[dim].push_back( x0 );
    while(x<=x_max)
    {

	derivs(x,Y,dY,dim);
	std::cout << x << " " << dY[0] << " " << dY[1] << std::endl;
	(*step)(Y, dY, dim, x, h,&(tmp_vec[0]), NULL, derivs, jacobn);
	ode_dat[dim].push_back( x += h );
	
	for(i=0; i<dim; i++)
	    ode_dat[i].push_back( Y[i] = tmp_vec[i] );
	
	ode_nstep++;
	h *= 1.002;
    }
}



void stiffint_adpt(vector<double> Y0, double x0, double x_max, double eps, double hi, int &nok, int &nbad,
	void (*derivs)(double, const vector<double> &, vector<double> &, int),
	void (*jacobn)(double, const vector<double> &, double**, int),
	void (*step)(vector<double>, vector<double>, int, double, double, vector<double> &, vector<double> &,
	    void(*)(double, const vector<double> &, vector<double> &, int),
	    void(*)(double, const vector<double> &, double**, int)),
	int order)
{
    const double safety = 0.9;
    const double pshrnk = -1.0/order;
    const double pgrow = -1.0/(order+1);
    const double errcon = pow(5.0/safety,1/pgrow);
    
    int i;
    double x=x0;
    double h=hi;
    double errmax;
    double htemp, hnext;
    double xnew;
    double hdid;
    int dim = Y0.size();

    vector<double> Y = Y0;
    vector<double> dY(dim);
    vector<double> tmp_vec(dim);
    vector<double> Yerr(dim);
    vector<double> Yscal(dim, 0.1);

    ode_nstep = 1;
    nok=nbad=0;
    
    ode_dat[dim].push_back( x0 );
    for(i=0; i<dim; i++)
	ode_dat[i].push_back( Y0[i] );

    while(x<=x_max)
    {

	derivs(x,Y,dY,dim);
	for(i=0; i<dim; i++)
	    Yscal[i]=fabs(Y[i])+fabs(dY[i]*h)+YMIN;

	hdid = h;
	for(;;)
	{
	    step(Y, dY, dim, x, h,tmp_vec, Yerr, derivs, jacobn);
	    errmax=0.0;
	    for(i=0;i<dim;i++)
		errmax = max(double(errmax),double(fabs(Yerr[i]/Yscal[i])));
	    errmax /= eps;
	    if(errmax <= 1.0) break;
	    htemp = safety*h*pow(errmax,pshrnk);
	    h = (h >= 0.0 ? max(htemp,0.1*h) : min(htemp,0.1*h));
	    xnew = x+h;
	    if (xnew == x) error("stepsize underflow");
	}
	if (errmax > errcon) hnext=safety*h*pow(errmax,pgrow);
	else hnext=5.0*h;//No more than a factor of 5 increase

	ode_dat[dim].push_back( x += h );
	for(i=0; i<dim; i++)
	    ode_dat[i].push_back( Y[i] = tmp_vec[i] );

	//x += h;
	if(h==hdid) nok++; else nbad++;
	h=hnext;
	if(ode_nstep++ > ode_maxlen)
	    exit(0);
    }
}



void odeint(double *Y0, int dim, double x0, double x_max, double eps, double hi, int *nok, int *nbad,
	void (*derivs)(double, double*, double*),
	void (*step)(double*, double*, int, double, double, double*, double*, 
	    void (*)(double, double*, double*)), int order)
{

    int i;
    
    double x=x0;
    double h;

    vector<double> Y(dim);
    vector<double> dY(dim);
    vector<double> tmp_vec(dim);

    for(i=0; i<dim; i++)
    {
	printf("%d\n",i);
	Y[i]=Y0[i];
    }


    h=hi;
    ode_nstep = 1;
    for(i=0; i<dim; i++)
	ode_dat[i][0] = Y0[i];
    ode_dat[dim][0] = x0;
    while(x<=x_max)
    {

	derivs(x,&(Y[0]),&(dY[0]));
	(*step)(&(Y[0]), &(dY[0]), dim, x, h,&(tmp_vec[0]), NULL, derivs);
	ode_dat[dim][ode_nstep] = x += h;
	
	for(i=0; i<dim; i++)
	    ode_dat[i][ode_nstep] = Y[i] = tmp_vec[i];
	
	ode_nstep++;

    }

}



void ode_init(int dim, int maxlen)
{
    ode_maxlen = maxlen;

    ode_tmp.resize(dim+1,2);
    
    ode_dat.resize(dim+1);
    for(int i=0; i<= dim; i++)
	ode_dat[i].reserve(maxlen);
	
    ode_jac.resize(dim, dim);
}


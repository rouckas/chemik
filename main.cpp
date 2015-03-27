#include "matrix.h"
#include "stiff.h"
#include "chemical.cpp"
#include <GetPot>

using namespace std;

double x0 = 0.0;
double x_max = 1e7;

double eps = 1e-3;

double hi = 1e-15;
int maxlen = 100000;

t_chemistry chemie;

void jacobn(double x, const vector<double> &y, double **dfdy, int n)
{
    for(int i=0; i<n; i++)
	chemie.set_concentration(i, y[i]);
    for(int i=0; i<n; i++)
	for(int j=0; j<n; j++)
	    dfdy[i][j] = chemie.jacobn(i,j);
}
void derivs(double x, const vector<double> &y, vector<double> &dydx, int n)
{
    for(int i=0; i<n; i++)
	chemie.set_concentration(i, y[i]);
    for(int i=0; i<n; i++)
	dydx[i] = chemie.deriv(i);
}


int main(int argc, char *argv[])
{

    GetPot cl(argc,argv);
    // load the configuration from file
    string config_file = cl("config","config.txt");
    GetPot config(config_file.c_str());

    string outfile = config("outfile", "res.dat");
    string specfile = config("specfile", "species.txt");
    string reactfile = config("reactfile", "reactions.txt");
    string plotfile = config("plotfile", "plotcmd.gnuplot");

    chemie.read_species(specfile.c_str());
    if( chemie.Nelem() == 0 )
    {
        cerr << "no elements loaded from species file '" << specfile << "', aborting" << endl;
        return 1;
    }
    chemie.read_reactions(reactfile.c_str());
    if( chemie.Nreactions() == 0 )
    {
        cerr << "no reactions loaded from reactions file '" << reactfile << "', aborting" << endl;
        return 1;
    }

    

    vector<double> Y0;
    for(int i=0; i<chemie.Nelem(); i++) Y0.push_back(chemie.Concentration(i));
    ode_init(chemie.Nelem(), 100000);


    int nok=0, nbad=0;
    stiffint_adpt(Y0, x0, x_max, eps, hi, nok, nbad, derivs, jacobn, rosenbrock_step, 4);
    cerr << nok << "  " << nbad << endl;


    chemie.sort_reactions();
    // Provnani celkoveho vlivu reakci
    
    
    for(int i=0; i<chemie.Nreactions(); i++)
    {
	//cout << setw(10) << chemie.Reaction(i)->importance5() << "  ";
	//chemie.Reaction(i)->print();
    }
    
    ofstream fw(outfile.c_str());
    fw << "#   ";
    for(int j=0; j<chemie.Nelem(); j++)
	fw << "  " << setw(7) << chemie.Species(j)->Name() ;
    fw << endl;
    for(int i=0;i<ode_nstep;i++)
    {
	fw << ode_dat[chemie.Nelem()][i];
	for(int j=0; j<chemie.Nelem(); j++)
	    fw << "  " << setw(7) << ode_dat[j][i] ;
	fw << endl;
    }
    fw.close();

    fw.open(plotfile.c_str());
    fw << "plot\\" << endl;
    for(int j=0; j<chemie.Nelem()-1; j++)
	fw << "  " << "'" << outfile << "' using 1:" << j+2 << " title '"
	    << chemie.Species(j)->Name() << "',\\" << endl ;
    fw << "  " << "'" << outfile << "' using 1:" << chemie.Nelem()+1 <<
	" title '" << chemie.Species(chemie.Nelem()-1)->Name() << "'" << endl ;
    fw.close();
    
    
    return 0;
}





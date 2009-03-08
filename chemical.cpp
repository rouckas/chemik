#ifndef CHEMICAL_H
#define CHEMICAL_H

#include "matrix.h"
#include <string>
#include <cmath>
#include <list>
#include <vector>
#include <algorithm>

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;


class t_reaction;

class t_species{
    private:
	std::string name;
	double concentration;
	vector<t_reaction*> is_product;
	vector<t_reaction*> is_reactant;
    public:
	t_species(std::string nam, double n0) 
	    : name(nam), concentration(n0) {}
	const string & Name() const
	    { return name; }
	double &Concentration() 
	    { return concentration; }
	void be_product(t_reaction * reaction)
	    { is_product.push_back(reaction); }
	void be_reactant(t_reaction * reaction)
	    { is_reactant.push_back(reaction); }
	double deriv();
	double jacobn(t_species *ps);
	void print()
	    { cout << name << " " << concentration << endl; }
	const vector<t_reaction*> & Is_product() const
	    { return is_product;};
	const vector<t_reaction*> & Is_reactant() const
	    { return is_reactant; };
	double rel_deriv(t_reaction * react);
	double rel_deriv_sum(t_reaction * react);
};

class t_reaction
{
    private:
	vector<t_species*> reactants;
	vector<t_species*> products;
	double rate_k;
    public:
	t_reaction(double Rate_k): rate_k(Rate_k) {};
	t_reaction(double Rate_k, const vector<t_species*> &Reactants, const vector<t_species*> &Products):
	    reactants(Reactants), products(Products), rate_k(Rate_k) {};
	double rate() { return rate_k; }
	void add_product(t_species * species)
	    { products.push_back(species); }
	void add_reactant(t_species * species)
	    { reactants.push_back(species); }
	double deriv()
	{
	    double res = rate();
	    for(unsigned int i=0; i<reactants.size(); i++)
		res *= reactants[i]->Concentration();
	    return res;
	}
	double jacobn(t_species *ps)
	{
	    double res = rate();
	    int order = 0;
	    for(unsigned int i=0; i<reactants.size(); i++)
	    {
		if( reactants[i] == ps )
		    if((order++)==0) continue;
		res *= reactants[i]->Concentration();
	    }
	    return res*order;
	}
	bool has_reactant(t_species *ps) const
	{
	    for(unsigned int i=0; i<reactants.size(); i++)
		if( ps == reactants[i] ) return true;
	    return false;
	}
	void print() const
	{
	    vector<t_species*>::const_iterator I = reactants.begin();
	    cout << rate_k << "\t\t";
	    if(reactants.size()>0)
		cout << (*I)->Name() ;
	    if(reactants.size()>1)
		for( I++; I != reactants.end(); I++)
		    cout << " " << (*I)->Name();

	    cout << " => ";

	    I = products.begin();
	    if(products.size()>0)
		cout << (*I)->Name() ;
	    if(products.size()>1)
		for( I++; I != products.end(); I++)
		    cout << " " << (*I)->Name();
	    cout << endl;
	}
	int nreactants(){ return reactants.size(); }
	double importance();
	double importance2();
	double importance3();
	double importance4();
	double importance5();
	static bool compare(t_reaction * a, t_reaction * b)
	{ return abs(a->importance5()) < abs(b->importance5()) ;}
};




class t_chemistry
{
    private:
	vector<t_species*> species;
	vector<t_reaction*> reactions;
    public:
	~t_chemistry()
	{
	    for(vector<t_species*>::iterator I=species.begin(); I != species.end(); I++)
		delete (*I);
	    for(vector<t_reaction*>::iterator I=reactions.begin(); I != reactions.end(); I++)
		delete (*I);
	}
	t_species * add_species(string nam, double n0)
	{
	    t_species *ps;
	    ps = new t_species(nam,n0);
	    species.push_back(ps);
	    return ps;
	}
	void set_concentration(int i, double concentration)
	{ species[i]->Concentration() = concentration; }
	t_reaction * add_reaction(double k, const vector<string> &reactants, const vector<string> &products)
	{
	    t_reaction *pr;
	    pr = new t_reaction(k);
	    reactions.push_back(pr);
	    vector<t_species*>::const_iterator I;
	    vector<string>::const_iterator J;
	    vector<t_species*> s_reactants, s_products;
	    //XXX dodelat kontrolu chyb
	    for(J=reactants.begin(); J !=reactants.end(); J++)
	    {
		for(I=species.begin(); I !=species.end(); I++)
		{
		    if( (*I)->Name() == *J )
		    {
			pr->add_reactant(*I);
			(*I)->be_reactant(pr);
			break;
		    }
		}
		if(I == species.end())
		{
		    cerr << "add_reaction: error: reactant " << (*J) << " not defined" << endl;
		}
	    }
	    for(J=products.begin(); J !=products.end(); J++)
	    {
		for(I=species.begin(); I !=species.end(); I++)
		{
		    if( (*I)->Name() == *J )
		    {
			pr->add_product(*I);
			(*I)->be_product(pr);
			break;
		    }
		}
	    }
		
	    return pr;
	}

	double jacobn(int i, int j)
	{ return species[i]->jacobn(species[j]); }
	
	double deriv(int i)
	{ return species[i]->deriv(); }
	
	double Concentration(int i)
	{ return species[i]->Concentration(); }
	
	int Nelem()
	{ return species.size(); }
	t_species* Species(int i)
	{ return species[i]; }

	int Nreactions()
	{ return reactions.size(); }
	t_reaction* Reaction(int i)
	{ return reactions[i]; }

	void read_reactions(const string fname)
	{


	    ifstream fr(fname.c_str());

	    string line;

	    vector<string> reactants;
	    vector<string> products;


	    double k;
	    string tmp;
	    t_reaction *pr;
	    while(fr.good())
	    {
		getline(fr, line);
		//cout << "good "<< fr.good() << " " << line << endl;

		istringstream s_line;
		s_line.str(line);

		if( !(s_line >> k) ) continue;
		while(s_line.good())
		{
		    s_line >> tmp;
		    if(tmp == "=>") break;
		    reactants.push_back(tmp);
		}
		while(s_line.good())
		{
		    s_line >> tmp;
		    products.push_back(tmp);
		}
		if( products.size()>0 || reactants.size()>0 )
		    pr = add_reaction(k, reactants, products);

		//pr->print();
		reactants.clear();
		products.clear();

	    }
	}

	void read_species(const string fname)
	{


	    ifstream fr(fname.c_str());

	    string line;

	    double conc;
	    string tmp;
	    t_species *ps;
	    while(fr.good())
	    {
		getline(fr, line);

		istringstream s_line;
		s_line.str(line);

		if( !(s_line >> tmp >> conc) || tmp[0]=='#')
		    continue;
		
		ps = add_species(tmp, conc);
		//ps->print();

	    }
	}
	void sort_reactions()
	{ sort(reactions.rbegin(), reactions.rend(), t_reaction::compare ); }
	
};



double t_species::deriv()
{
    double res=0;
    
    for(unsigned int i=0; i<is_product.size(); i++)
	res += is_product[i]->deriv();
    for(unsigned int i=0; i<is_reactant.size(); i++)
	res -= is_reactant[i]->deriv();
    return res;
}
double t_species::jacobn(t_species *ps)
{
    double res=0;
    unsigned int i;
    for(i=0; i<is_product.size(); i++)
	if( is_product[i]->has_reactant(ps) )
	    res += is_product[i]->jacobn(ps);
    for(i=0; i<is_reactant.size(); i++)
	if( is_reactant[i]->has_reactant(ps) )
	    res -= is_reactant[i]->jacobn(ps);
    return res;
}
double t_species::rel_deriv(t_reaction * react)
{
    double res = 0;
    double der = 0;
    double derr = react->deriv();// derivace prislusejici zkoumane rovnici
    double tmp = 0;
    //spocteme derivaci s kladnym znaminkem
    for(unsigned int i=0; i<is_product.size(); i++)
    {
	der += is_product[i]->deriv();
	if(is_product[i] == react)
	    tmp += derr;
    }
    res = tmp/der;
    der = 0;
    tmp = 0;
    //spocteme derivaci se zapornym znaminkem
    for(unsigned int i=0; i<is_reactant.size(); i++)
    {
	der -= is_reactant[i]->deriv();
	if(is_reactant[i] == react)
	    tmp -= derr;
    }
    return res+tmp/der;
}
double t_species::rel_deriv_sum(t_reaction * react)
{
    double res = 0;
    double der = 0;
    double derr = react->deriv();// derivace prislusejici zkoumane rovnici
    double tmp = 0;
    //spocteme derivaci s kladnym znaminkem
    for(unsigned int i=0; i<is_product.size(); i++)
    {
	der += is_product[i]->deriv();
	if(is_product[i] == react)
	    tmp += derr;
    }
    res = tmp/der;
    der = 0;
    tmp = 0;
    //spocteme derivaci se zapornym znaminkem
    for(unsigned int i=0; i<is_reactant.size(); i++)
    {
	der -= is_reactant[i]->deriv();
	if(is_reactant[i] == react)
	    tmp -= derr;
    }
    return res-tmp/der;
}





double t_reaction::importance()
{
    double res = 0;
    double der = deriv();
    for(unsigned int i=0; i<reactants.size(); i++)
	res += der / reactants[i]->Concentration();
    for(unsigned int i=0; i<products.size(); i++)
	res += der / products[i]->Concentration();
    return res;
}
double t_reaction::importance2()
{
    double res = 0;
    for(unsigned int i=0; i<reactants.size(); i++)
    {
	unsigned int j;
	for(j=0; j<i; j++)
	    if( reactants[i] == reactants[j] )
		break;
	if(j==i)	//pricitame jen unikatni slozky
	    res += reactants[i]->rel_deriv(this);
    }
    for(unsigned int i=0; i<products.size(); i++)
    {
	unsigned int j;
	for(j=0; j<i; j++)
	    if( products[i] == products[j] )
		break;
	if(j==i)
	    res += products[i]->rel_deriv(this);
    }
    return res;
}
double t_reaction::importance3()
{
    double res = 0;
    double der = deriv();
    for(unsigned int i=0; i<reactants.size(); i++)
	res += der / reactants[i]->Concentration();
    for(unsigned int i=0; i<products.size(); i++)
	res -= der / products[i]->Concentration();
    return res;
}
double t_reaction::importance4()
{
    double res = 0;
    for(unsigned int i=0; i<reactants.size(); i++)
    {
	unsigned int j;
	for(j=0; j<i; j++)
	    if( reactants[i] == reactants[j] )
		break;
	if(j==i)	//pricitame jen unikatni slozky
	    res += reactants[i]->rel_deriv(this);
    }
    for(unsigned int i=0; i<products.size(); i++)
    {
	unsigned int j1, j2;
	for(j1=0; j1<reactants.size(); j1++)
	    if( products[i] == reactants[j1] )
		break;
	for(j2=0; j2<i; j2++)
	    if( products[i] == products[j2] )
		break;
	if(j2==i && j1==reactants.size())
	    res += products[i]->rel_deriv(this);
    }
    return res;
}
double t_reaction::importance5()
{
    double res = 0;
    for(unsigned int i=0; i<reactants.size(); i++)
    {
	unsigned int j;
	for(j=0; j<i; j++)
	    if( reactants[i] == reactants[j] )
		break;
	if(j==i)	//pricitame jen unikatni slozky
	    res += reactants[i]->rel_deriv_sum(this);
    }
    for(unsigned int i=0; i<products.size(); i++)
    {
	unsigned int j1, j2;
	for(j1=0; j1<reactants.size(); j1++)
	    if( products[i] == reactants[j1] )
		break;
	for(j2=0; j2<i; j2++)
	    if( products[i] == products[j2] )
		break;
	if(j2==i && j1==reactants.size())
	    res += products[i]->rel_deriv_sum(this);
    }
    return res;
}
#endif

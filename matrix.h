
#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <vector>
using namespace std;


template<class T> class t_matrix
{
    protected:
	T *m;
	T **vm;
	int cols, rows;

    public:
	t_matrix() : cols(0), rows(0) {} ;
	t_matrix(int r, int c);

	~t_matrix()
	{ delete [] m;  delete [] vm; }

	T * operator [] (int i) const
	{ return vm[i];	}

	int Rows() const { return rows; }
	int Cols() const { return cols; }

	t_matrix<T> & operator = (const t_matrix<T> &r);
	void resize(int r, int c);

	t_matrix (const t_matrix<T> &r);

	operator T** ()
	{ return vm; }
};


template<class T> class t_numatrix : public t_matrix<T>
{
    protected:
	t_matrix<T>::rows;
	t_matrix<T>::cols;
	t_matrix<T>::vm;
	t_matrix<T>::m;
    public:
	t_numatrix() : t_matrix<T>() {} ;
	t_numatrix(int rows, int cols) : t_matrix<T>(rows, cols) {} ;
	t_numatrix(t_numatrix<T> &r) : t_matrix<T>(r) {} ;
    
    t_numatrix<T> operator * (const t_numatrix &r);

    void print();

};
	
/*
template<class T> t_vector<T> operator * (const t_matrix<T> &l, const t_vector<T> &r);

template<class T> t_vector<T> operator * (const t_vector<T> &l, const t_matrix<T> &r);
*/
template<class T> class t_LUmatrix : public t_numatrix<T>
{
    protected:
	t_matrix<T>::rows;
	t_matrix<T>::cols;
	t_matrix<T>::vm;
	t_matrix<T>::m;
    private:
	vector<int> indx;
	signed char sign;

    public:
	t_LUmatrix() : t_numatrix<T>() {} ;
	t_LUmatrix(int rows, int cols);

	t_LUmatrix<T> & operator = ( const t_numatrix<T> &r );
	void resize(int r, int c);

	void decompose ();

	void LUbksb(vector<T> &b);
};

template class t_matrix<int>;

template class t_matrix<double>;
template class t_numatrix<double>;
template class t_LUmatrix<double>;



static void error(const char *str)
{
    fprintf(stderr, "%s\n", str);
    exit(1);
}


template<class T>
t_matrix<T>::t_matrix(int r, int c)
{
    cols = c;
    rows = r;

    m = new T [cols*rows];
    vm = new T* [rows];

    T *p1;
    p1 = m;
    for(int i=0; i<cols*rows; i++)
	*(p1++) = T();		//inicializace prvku

    for(int i=0; i<rows; i++)
	vm[i] = &m[i*cols];
}

template<class T>
t_matrix<T> & t_matrix<T>::operator = (const t_matrix<T> &r)
{
    if (cols != r.cols || rows != r.rows)
	error("matrix: operator = : size mismatch error");

    for(int i=0; i<cols*rows; i++)
	m[i] = r.m[i];

    return *this;
}

template<class T>
void t_matrix<T>::resize(int r, int c)
{
    cols = c;
    rows = r;

    delete [] m;
    delete [] vm;
    m = new T [cols*rows];
    vm = new T* [rows];

    T *p1;
    p1 = m;
    for(int i=0; i<cols*rows; i++)
	*(p1++) = T();		//inicializace prvku

    for(int i=0; i<rows; i++)
	vm[i] = &m[i*cols];
}

template<class T>
t_matrix<T>::t_matrix (const t_matrix<T> &r)
{
    cols = r.cols;
    rows = r.rows;

    m = new T [cols*rows];
    vm = new T* [rows];

    for(int i=0; i<cols*rows; i++)
	m[i] = r.m[i];

    for(int i=0; i<rows; i++)
	vm[i] = &m[i*cols];
}




template<class T>
t_numatrix<T> t_numatrix<T>::operator * (const t_numatrix &r)
{
    t_numatrix<T> &l = *this;
    if(l.cols != r.rows)
	error("t_numatrix: operator * : incompatible matrices");

    int row = l.rows;
    int col = r.cols;
    t_numatrix<T> prod(row, col);

    T sum;
    int i, j, k;
    for(i=0; i<col; i++)
	for(j=0; j<row; j++)
	{
	    sum = T(0);
	    for(k=0; k<l.rows; k++)
		sum += l[i][k]*r[k][j];
	    prod[i][j] = sum;
	}

    return prod;
}

template<class T>
void t_numatrix<T>::print()
{
    for(int i=0; i<rows; i++)
    {
	for(int j=0; j<cols; j++)
	    std::cout << std::setw(5) << vm[i][j] << "  ";
	std::cout << std::endl;
    }
}

	
/*
template<class T> t_vector<T> operator * (const t_matrix<T> &l, const t_vector<T> &r)
{
    int dim = l.Rows();
    t_vector<T> res(dim);
    int i, j;
    T sum;
    for(i=0; i<dim; i++)
    {
	sum = T(0);
	for(j=0; j<l.Cols(); j++)
	    sum += l[i][j]*r[j];
	res[i] = sum;
    }
    return res;
}

template<class T> t_vector<T> operator * (const t_vector<T> &l, const t_matrix<T> &r)
{
    int dim = r.Cols();
    t_vector<T> res(dim);
    int i, j;
    T sum;
    for(i=0; i<dim; i++)
    {
	sum = T(0);
	for(j=0; j<r.Rows(); j++)
	    sum += l[j]*r[j][i];
	res[i] = sum;
    }
    return res;
}
*/
template<class T> T abs(const T x)
{
    return x>T(0) ? x : -x;
}

template<class T>
t_LUmatrix<T>::t_LUmatrix(int rows, int cols) : t_numatrix<T>(rows, cols), indx(rows)
{
    if(rows != cols)
	error("t_LUmatrix: matrix not rectangular");
}

template<class T>
t_LUmatrix<T> & t_LUmatrix<T>::operator = ( const t_numatrix<T> &r )
{
    if(rows != r.Rows() || rows != r.Cols())
	error("t_LUmatrix: operator = : size mismatch");
    this->t_numatrix<T>::operator=(r);
    return *this;
}
	    
template<class T>
void t_LUmatrix<T>::resize(int r, int c)
{
    if(r != c)
	error("t_LUmatrix: matrix not rectangular");
    cols = c;
    rows = r;

    indx.resize(r);
    
    delete [] m;
    delete [] vm;
    m = new T [cols*rows];
    vm = new T* [rows];

    T *p1;
    p1 = m;
    for(int i=0; i<cols*rows; i++)
	*(p1++) = T();		//inicializace prvku

    for(int i=0; i<rows; i++)
	vm[i] = &m[i*cols];
}


template<class T>
void t_LUmatrix<T>::decompose ()
{
#define TINY 1.0e-100
    t_LUmatrix<T> &a = *this;

    int i, imax=0, j, k;
    int n = cols;
    T big, dum, sum, temp;

    vector<T> vv(cols);

    sign = 1;
    for(i=0; i<n; i++)
    {
	big = T(0);
	for(j=0; j<n; j++)
	    if((temp=abs<T>(a[i][j])) > big) 
		big = temp;
	if(big==T(0))
	    error("ludcmp: singular matrix");
	vv[i] = T(1)/big;
    }

    for(j=0; j<n; j++)
    {
	for(i=0; i<j; i++)
	{
	    sum = a[i][j];
	    for(k=0; k<i; k++)
		sum -= a[i][k]*a[k][j];
	    a[i][j] = sum;
	}
	big = T(0);
	for(i=j; i<n; i++)
	{
	    sum = a[i][j];
	    for(k=0; k<j; k++)
		sum -= a[i][k]*a[k][j];
	    a[i][j] = sum;
	    if( (dum=vv[i]*abs<T>(sum)) >= big )
	    {
		big = dum;
		imax = i;
	    }
	}
	if(j != imax)
	{
	    for(k=0; k<n; k++)
	    {
		dum = a[imax][k];
		a[imax][k] = a[j][k];
		a[j][k] = dum;
	    }
	    sign = -sign;
	    vv[imax] = vv[j];
	}
	indx[j] = imax;
	if(a[j][j] == T(0))
	    a[j][j] = TINY;
	if(j != n)
	{
	    dum = T(1)/(a[j][j]);
	    for(i=j+1; i<n; i++)
		a[i][j] *= dum;
	}
    }
}

template<class T>
void t_LUmatrix<T>::LUbksb(vector<T> &b)
{
    int i, ii=-1, ip, j;
    int n = cols;
    t_LUmatrix<T> &a = *this;
    T sum;

    for(i=0; i<n; i++)
    {
	ip = indx[i];
	sum = b[ip];
	b[ip] = b[i];
	if(ii>=0)
	    for(j=ii; j<i; j++)
		sum -= a[i][j]*b[j];
	else if(sum != T(0))
	    ii = i;
	b[i] = sum;
    }
    for(i=n-1; i>=0; i--)
    {
	sum = b[i];
	for(j=i+1; j<n; j++)
	    sum -= a[i][j]*b[j];
	b[i] = sum/a[i][i];
    }
}
#endif

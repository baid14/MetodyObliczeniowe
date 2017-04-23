#include "stdafx.h"
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;

int N = 10;
const double x1 = 0.0;
const double x2 = 2.0;

double Conventional(double *L, double *D, double *U, double *B, double h);
double Numerical(double *L, double *D, double *U,double *B, double h);
double Analitical(double x);
double Error(double *blad);
double Precision(double x1, double y1, double x2, double y2);
double *Thomas(double *l, double *d, double *u, double *B);
double *Transform_matrix(double *l, double *d, double *u);
double *Solve_Thomas(double *l, double *ni, double *u, double *B);

int _tmain(int argc, _TCHAR* argv[]) {
	ofstream f;
	f.open("wyniki.txt");
	f.setf( ios::scientific );
	
	cout<< "Rzedy Dokladnosci\nMetoda Konwencjonalna = " << Precision(2.222222e-001, 5.421619e-005, 1.005025e-002, 1.116594e-007) <<"\nMetoda Numerowa = " << Precision(2.222222e-001, 3.350231e-008, 2.898551e-002, 9.752199e-012) <<"\n" ;

	double *L = new double[20000];
	double *D = new double[20000];
	double *U = new double[20000];
	double *B = new double[20000];

	for(int i=0 ; i<2000 ; i++) {
		double h = abs(x2-x1)/(N-1);
		f << h << "\t" << Conventional(L,D,U,B,h) << "\t" << Numerical(L,D,U,B,h) << endl;
		N += 10;
	}

	delete []L;
	delete []D;
	delete []U;
	delete []B;

	system("Pause");
	return 0;
}
double Conventional(double *L, double *D, double *U,double *B, double h) {
	double *error = new double[N];
	double *x = new double[N];
	double a = 0.0;

	U[0] = 0.0;
	D[0] = 1.0;
	L[0] = 0.0;
	B[0] = 1.0;

	for(int i=1 ; i<N-1 ; i++) {
		L[i] = 1.0/(h*h);
		D[i] = -(8.0+h*h)/(4.0*h*h);
		U[i] = 1.0/(h*h);
		B[i] = 0.0;
	}

	U[N-1] = 0.0;
	D[N-1] = 1.0;
	L[N-1] = 0.0;
	B[N-1] = 0.0;

	x = Thomas(L,D,U,B);
	for(int i=0; i<N; i++) {
		error[i] = abs(x[i] - Analitical(a));
		a += h;		
	}
	delete[] x;
	return Error(error);
}
double Numerical(double *L, double *D, double *U,double *B, double h) {
	double *error = new double[N];
	double *x =new double[N];
	double a = 0.0;
	
	U[0] = 0.0;
	D[0] = 1.0;
	L[0] = 0.0;
	B[0] = 1.0;
	
	for(int i = 1; i < N-1; i++) {
		L[i] = (48.0-h*h)/(48.0*h*h);
		D[i] = (-48.0-5.0*h*h)/(24.0*h*h);
		U[i] = (48.0-h*h)/(48.0*h*h);
		B[i] = 0.0;
	}

	U[N-1] = 0.0;
	D[N-1] = 1.0;
	L[N-1] = 0.0;
	B[N-1] = 0.0;

    x = Thomas(L,D,U,B);

	for(int i=0; i<N; i++) {
		error[i] = abs(x[i] - Analitical(a));
		a+=h;
	}
	delete[] x;
	return Error(error);
}
double Analitical(double x) {
	return (exp(x/2.0) - exp(2.0-x/2.0))/(1.0-exp(2.0));
}
double Error(double *error) {
	double max = 0.0;
    for(int k=0 ; k<N ; k++)
		if(abs(error[k])> max)
			max = abs(error[k]);
	delete[] error;
	return max;
}
double Precision(double x1, double y1, double x2, double y2) {
	return ( (log(y1)/log(10.0) - log(y2)/log(10.0) ) / (log(x1)/log(10.0) - log(x2)/log(10.0)) );
}
double *Thomas(double *l, double *d, double *u, double *B) {
	return Solve_Thomas(l, Transform_matrix(l,d,u), u, B);
}
double *Transform_matrix(double *l, double *d, double *u) {
	double *ni = new double[N];
	ni[0] = d[0];
	for(int i = 1; i<N; i++)
		ni[i] = d[i] - (l[i]*u[i-1])/ni[i-1];
	return ni;
}
double *Solve_Thomas(double *l, double *ni, double *u, double *B) {
	double *b = new double[N];
	double *r = new double[N];
	r[0] = B[0];
	for(int i=1 ; i<N ; i++)
		r[i] = B[i] - (l[i]*r[i-1])/ni[i-1];
	b[N-1] = r[N-1]/ni[N-1];
	for(int i=N-2 ; i>=0 ; i--)
		b[i] = (r[i] - u[i]*b[i+1])/ni[i];
	delete[] r;
	delete[] ni;
	return b;
}
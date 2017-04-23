#pragma once
#include "stdafx.h"

class ProjectD
{
public:
	ProjectD(void);
	~ProjectD(void);

	void info();
	void compute(int n);
	void setMethod(int meth);
	long double exactSolution(int t, int x)
	{
		return 1.0 - (r/pos[x]) * erfc((pos[x]-r)/(2.0*sqrt(D*times[t])));
	}
	bool computed()
	{
		if ( N > 0 ) 
			return true;
		else
			return false;
	}
	int getN(){	return N;}
	long double getH(){ return h;}
	int getM(){	return M;}
	long double getDt(){ return dt;}
	long double getA(){ return a;}
	long double getR(){ return r;}
	long double** getX(){ return x;}
	long double *getTimes() {return times;}
	long double *getPos() {return pos;}
	long double getElapsed(){ return elapsed;}
	std::string getName(){ return methName;}

private:
	int D;
	long double TMAX;
	int N;				//ilosc podzialow przedzialu przestrzennego (zmienna x)
	long double h;			//krok przestrzenny
	int M;				//ilosc podzialow przedzialu czasowego (zmienna t)
	long double dt;			//krok czasowy;
	long double a;
	long double r;
	long double lambda;		
	long double lambda2;
	long double *pos;
	long double *times;

	long double elapsed;

	std::string methName;

	long double **x;			//macierz rozwiazan zagadnienia (serii ukladow rownan Ax = b)

	void computeCNT(int n);
	void computeKMB(int n);

	void initialValues(long double *b);
	
	void diskretCNT(long double *l, long double *d, long double *u); 
	void diskretLT(long double *l, long double *d, long double *u); 

	void thomasReduce(long double *l, long double *d, long double *u);        
	void thomasSolve(long double *l, long double *d, long double *u, long double *x, long double *b);

	void reset();
};
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
	double exactSolution(int t, int x)
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
	double getH(){ return h;}
	int getM(){	return M;}
	double getDt(){ return dt;}
	double getA(){ return a;}
	double getR(){ return r;}
	double** getX(){ return x;}
	double *getTimes() {return times;}
	double *getPos() {return pos;}
	double getElapsed(){ return elapsed;}
	std::string getName(){ return methName;}

private:
	int D;
	double TMAX;
	int N;				//ilosc podzialow przedzialu przestrzennego (zmienna x)
	double h;			//krok przestrzenny
	int M;				//ilosc podzialow przedzialu czasowego (zmienna t)
	double dt;			//krok czasowy;
	double a;
	double r;
	double lambda;		
	double lambda2;
	double *pos;
	double *times;

	double elapsed;

	std::string methName;

	double **x;			//macierz rozwiazan zagadnienia (serii ukladow rownan Ax = b)

	void computeCNT(int n);
	void computeKMB(int n);

	void initialValues(double *b);
	
	void diskretCNT(double *l, double *d, double *u); 
	void diskretLT(double *l, double *d, double *u); 

	void thomasReduce(double *l, double *d, double *u);        
	void thomasSolve(double *l, double *d, double *u, double *x, double *b);

	void reset();
};
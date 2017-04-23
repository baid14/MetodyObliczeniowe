#include <iostream>
#include <cmath>
#include "calerf.h"
#include <fstream>
using namespace std;

const double TMIN=0;
const double TMAX=2;
const double D=1;
const double h=0.1;
const double LAMBDA=0.4;
const double dt=(LAMBDA*(h*h))/D;
const double a=6*sqrt(D*TMAX);
const int n=(int)(TMAX/dt)+1;
const int m=(int)a/h;

void analit(double **Analit, int n , int m, fstream & file)
{
     double x=h;
     double t=dt;
     for(int j=0;j<m;j++)
             Analit[0][j]=0;
             
     for(int i=0;i<n;i++)
             Analit[i][m-1]=0;
             
     for(int i=0;i<n;i++)
             Analit[i][0]=1;
             
     for(int i=1;i<n;i++)
     {
             for(int j=1;j<m-1;j++)
             {
                     Analit[i][j]=erfc(x/(2*sqrt(D*t)));
                     x=x+h;
             }
             x=h;
             t=t+dt;
     }
     
     for(int i=0;i<n;i++)
     {
             for(int j=0;j<m;j++)
             {
                     file<<Analit[i][j]<<";";
             }
             file<<"\n";
     }
     
}

void bezpos(double **A, int n, int m, fstream &file2)
{
     for(int j=0;j<m;j++)
             A[0][j]=0;
             
     for(int i=0;i<n;i++)
             A[i][m-1]=0;
             
     for(int i=0;i<n;i++)
             A[i][0]=1;
             
     for(int i=1;i<n;i++)
     {
             for(int j=1;j<m-1;j++)
             {
                     A[i][j]=A[i-1][j]+LAMBDA*(A[i-1][j-1]-(2*A[i-1][j])+A[i-1][j+1]);             
             }
     }
     
     for(int i=0;i<n;i++)
     {
             for(int j=0;j<m;j++)
             {
                     file2<<A[i][j]<<";";
             }
             file2<<"\n";
     }
}

void blad(double **Blad, double **Analit, double **A, int n, int m, fstream &file3)
{
     for(int i=0;i<n;i++)
     {
             for(int j=0;j<m;j++)
             {
                     Blad[i][j]=fabs(Analit[i][j]-A[i][j]);
             }
     }
     
     for(int i=0;i<n;i++)
     {
             for(int j=0;j<m;j++)
             {
                     file3<<Blad[i][j]<<";";
             }
             file3<<"\n";
     }
                     
} 
int main()
{
    fstream file,file2,file3;
    file.open("Analit.csv",ios::out);
    file2.open("Numer.csv",ios::out);
    file3.open("Blad.csv",ios::out);
    double **Analit;
    double **A;
    double **Blad;
    
    Analit =new double* [n];
	for(int i=0;i<n;i++)
		Analit[i]=new double[m];
		
	A =new double* [n];
	for(int i=0;i<n;i++)
		A[i]=new double[m];
			
	Blad =new double* [n];
	for(int i=0;i<n;i++)
		Blad[i]=new double[m];
		
		analit(Analit,n,m,file);
 	    bezpos(A,n,m,file2);
 	    blad(Blad,Analit,A,n,m,file3);
 	    
	file.close();
    file2.close();
    file3.close();
    
	for(int i=0; i<n; i++)
    {
          delete [] Analit[i];
          delete [] A[i];
          delete [] Blad[i];
    }
    delete [] Analit;
    delete [] A;
    delete [] Blad;
	Analit=A=Blad=NULL;	
    system("pause");
    return 0;
}

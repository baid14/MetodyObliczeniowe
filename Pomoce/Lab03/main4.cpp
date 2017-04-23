#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<windows.h>

/*

Metody Obliczeniowe
By Wahuu

*/


double eps = 0.0000001;

double **nowa_macierz(int, int); 
void usun_macierz(double **, int); 
void funkcje(double *, double *);
void macierz_jacobiego(double *, double **);
void obl_poprawki(double **, double *, double *); 
void obl_n_pierw(double *, double *);
double norma_macierzy(double *);


int main() 

{
	int i = 0;
	double **macierz = nowa_macierz(2,2);
	double *pierw = new double[2];
	double *popr = new double[2];
	double *fun = new double[2];

	//Podanie danych z zadania
	pierw[0] = 0.0;
	pierw[1] = sqrt(8.0);

	printf("METODY OBLICZENIOWE Zad 4\n\n");
	printf("****************************************************\n\n");
    printf("+-----+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+\n");
    printf("|  i  |        x        |        y        |     F1(x,y)     |      F2(x,y)    |   Poprawka 1    |    Poprawka 2   |\n");
    printf("+-----+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+\n");
    
	while(1) 
    {

		funkcje(fun,pierw);
		if(norma_macierzy(fun) < eps) 
              break;
              
  		macierz_jacobiego(pierw, macierz);
		obl_poprawki(macierz,popr,fun);

		obl_n_pierw(pierw, popr);
		i++; 
		if(i > 53 || norma_macierzy(popr) < eps) 
              break;  
			      	
        printf("| %2d. | % 15.11lf | % 15.11lf | % 15.11lf | % 15.11lf | % 15.11lf | % 15.11lf |\n", i, pierw[0], pierw[1], fun[0], fun[1],popr[0],popr[1]);

	}
    printf("+-----+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+\n");
	printf("\n\nRozwiazanie:\n   x1=%e\n   x2= %e\n", pierw[0], pierw[1]);
	usun_macierz(macierz, 2);
    getchar();
    return 0;

}



double **nowa_macierz(int n, int m) 
{                       
	int i;
	double **a;
	a = new double*[n];
	for(i=0; i<n; i++) 
         a[i] = new double[m];
	return a;
}



void usun_macierz(double **a, int n) 
{                   
	int i;
	for(i=n-1; i>=0; i--) 
         delete []a[i];

	delete []a;
}



void funkcje(double *funkcje, double *pierw) 
{                       
	funkcje[0] = pierw[1]*pierw[1]*pierw[1] + 4.0f*pierw[0];
	funkcje[1] = pierw[0]*pierw[0] + pierw[1]*pierw[1] - 8.0;
}



void macierz_jacobiego(double *pierw, double **jac) 
{                       
	jac[0][0] = 4.0;
	jac[0][1] = 3.0*pierw[1]*pierw[1];
	jac[1][0] = 2.0*pierw[0];
	jac[1][1] = 2.0*pierw[1];
}



void obl_poprawki(double **jacobi, double *popr, double *fun) 
{			
	popr[0] = (fun[0]*jacobi[1][1] - jacobi[0][1]*fun[1]) / (jacobi[0][0]*jacobi[1][1] - jacobi[0][1]*jacobi[1][0]);
	popr[1] = (jacobi[0][0]*fun[1] - fun[0]*jacobi[1][0]) / (jacobi[0][0]*jacobi[1][1] - jacobi[0][1]*jacobi[1][0]);
}



void obl_n_pierw(double *pierw, double *popr) 
{                    
	pierw[0] = pierw[0] - popr[0];
	pierw[1] = pierw[1] - popr[1];
}



double norma_macierzy(double *popr) 
{                                   
	if(fabs(popr[0]) < fabs(popr[1])) 
          return fabs(popr[1]);
	else 
          return fabs(popr[0]);                               
}





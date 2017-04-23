#include <iostream>
#include <math.h>
//#include <fstream>
#include <Windows.h>

/*

Metody Obliczeniowe
By Wahuu

*/



using namespace std;

double EPS = 0.000001;

double **create_matrix(int n, int m);
void delete_matrix(double **a, int n);
void fill_matrix(double **a, int n);
void print_vector(double *b,int n);
double vect_norm(double *v);
double* residuum(double **a, double *x, double *b, int n);
void jacobi(double **matrix, double *b, double *x,  int n);
void gauss_seidel(double **matrix, double *b, double *x, int n);
void sor(double **matrix, double *b, double *x, double omega, int n);

double **create_matrix(int n, int m)
{
	double **a = new double*[n];
	for(int i=0; i<n; i++)
	    a[i] = new double[m];
	return a;
}

void delete_matrix(double **a, int n)
{
	for(int i=n-1; i>=0; i--)
        delete []a[i];
	delete []a;
}

void fill_matrix(double **A)
{
	A[0][0] = 10.0; A[0][1] = 1.0; A[0][2] = 2.0;
	A[1][0] = 1.0; A[1][1] = 20.0; A[1][2] = 5.0;
	A[2][0] = 3.0; A[2][1] = 4.0; A[2][2] = 30.0; ;
}

void print_vector(double *b,int n)
{
     for(int i=0; i<n; i++)
     {
              cout << b[i] << endl;
     }
}
/*
double vect_norm(double *v)
{
	if(fabs(v[0]) < fabs(v[1]))
        return fabs(v[1]);
	else
      return fabs(v[0]);
}
//r=Ax-b
double* residuum(double **a, double *x, double *b, int n)
{
     double sum = 0;
     double *r = new double[n];

     for(int i=0; i<n; i++)
     {
              for(int j=0; j<n; j++)
              {
                       sum = sum + a[i][j] * x[j];
              }
              r[i] = sum - b[i];
              sum = 0;
     }
     return r;
}
*/

double est( double *x, double *x_nowe )
{

	x[0] = x[0] - x_nowe[0];
	x[1] = x[1] - x_nowe[1];
	x[2] = x[2] - x_nowe[2];

	double norma = x[0];
	if ( x[1] > norma ) norma = x[1];
	if ( x[2] > norma ) norma = x[2];


	return norma;
}

double redisuum( double **A, double *b, double *x_nowe )
{
	double Ax[3];

	Ax[0] = fabs( (A[0][0] * x_nowe[0] + A[0][1] * x_nowe[1] + A[0][2] * x_nowe[2]) - b[0] );
	Ax[1] = fabs( (A[1][0] * x_nowe[0] + A[1][1] * x_nowe[1] + A[1][2] * x_nowe[2]) - b[1] );
	Ax[2] = fabs( (A[2][0] * x_nowe[0] + A[2][1] * x_nowe[1] + A[2][2] * x_nowe[2]) - b[2] );

	double norma = Ax[0];
	if ( Ax[1] > norma ) norma = Ax[1];
	if ( Ax[2] > norma ) norma = Ax[2];


	return norma;
}

void jacobi(double **matrix, double *b, double *x,  int n)
{
	int iter = 0;
	double *x1 = new double[n];
	double *Xn;                                         // wektor przyblizen
	double *Xp;                                       // przyblizenie poczatkowe
	double s = 0;

 	while(1)
	{
           if(iter == (iter / 2) * 2)
           {
                                                    Xn = x1;
                                                    Xp = x;
           }
           else
           {
                     Xn = x;
                     Xp = x1;
           }

           for(int i=0; i<n; i++)
           {
                   for(int j=0; j<n; j++)
                   {
                            if(j != i)
                                s = s + matrix[i][j] * Xp[j];
                   }
                   Xn[i] = (1.0 / matrix[i][i]) * (b[i] - s);
                   s = 0;
           }
         //  norma wartosci bezwzgl f. dla kolejnego przybliz || roznica miedzy kolenymi przyblizeniami || skonczona ilosc iteracji
           if(redisuum(matrix, b, Xn)  < EPS || est(Xp, Xn) < EPS || iter > 53)
               break;
           iter++;
		   cout <<"i: "<< iter << endl;
		   print_vector(Xn, 3);
           cout<<"Residuum: " << redisuum(matrix, b, Xn) << endl << endl;
    }
}

void gauss_seidel(double **matrix, double *b, double *x, int n)
{
	int iter = 1;
	double *Xp = new double[n];
	double s = 0;

    while(1)
	{
		for(int i=0; i<n; i++)
		{
			Xp[i] = x[i];
            //L
			for(int j=0; j<=i-1; j++)
			{
				s = s + matrix[i][j] * x[j];
			}
			//U
			for(int j=i+1; j<n; j++)
			{
				s = s + matrix[i][j] * x[j];
			}

			x[i] = (1.0 / matrix[i][i]) * (b[i] - s);
			s = 0;
		}
        //norma wartosci bezwzgl f. wekt resudialnego dla kolejn przybliz || roznica miedzy kolenymi przyblizeniami || skonczona ilosc iteracji
	    if(fabs((vect_norm(residuum(matrix, x, b, 2)))) < EPS || fabs((vect_norm(x) - vect_norm(Xp))) < EPS || iter > 53)
               break;
        cout <<"i: "<< iter << endl;
        print_vector(x, 2);
        cout<<"Residuum: " << fabs(vect_norm(residuum(matrix, x, b, 2))) << endl << endl;;
        iter++;
	}
}

void sor(double **matrix, double *b, double *x, double omega, int n)
{
    int iter = 1;
	double *Xp = new double[n];
	double s = 0;

    while(1)
	{
		for(int i=0; i<n; i++)
		{
			Xp[i] = x[i];
			//L
			for(int j=0; j<=i-1; j++)
			{
				s = s + matrix[i][j] * x[j];
			}
			//D/U
			for(int j=i; j<n; j++)
			{
                       if(j == i)
                           s = s + (1 - (1.0 / omega)) * matrix[i][i] * x[j];
                       else
                          s = s + matrix[i][j] * x[j];
			}

			x[i] = (1.0 / matrix[i][i]) * (b[i] - s) * omega;
			s = 0;
		}
        //norma wartosci bezwzgl f. wekt resudialnego dla kolejn przybliz || roznica miedzy kolenymi przyblizeniami || skonczona ilosc iteracji
	    if(fabs((vect_norm(residuum(matrix, x, b, 2)))) < EPS || fabs((vect_norm(x) - vect_norm(Xp))) < EPS || iter > 53)
               break;
        cout <<"i: "<< iter << endl;
        print_vector(x, 2);
        cout<<"Residuum: " << fabs(vect_norm(residuum(matrix, x, b, 2))) << endl << endl;
        iter++;
	}
}

int main()
{
      double **matrix = create_matrix(3, 3);
      double b[3];
      b[0] = 8.0; b[1] = -4.0; b[2] = -27.0;        // wektor b
      double x[3];
      x[0] = 2.0; x[1] = 2.0; x[2] = 2.0;         // przyblizenie poczatkowe
      int c;

      fill_matrix(matrix);
      cout << "Wybierz metode:" << endl;
      cout << "1. Jacobi" << endl;
      cout << "2. Gauss-Seidel"<<endl;
      cout << "3. SOR" << endl;
      cin >> c;

      switch(c)
      {
               case 1:
                    system("CLS");
                    jacobi(matrix, b, x, 3);
                    system("pause");
                    break;
              /* case 2:
                    system("CLS");
                    gauss_seidel(matrix, b, x, 2);
                    system("pause");
                    break;
               case 3:
                    system("CLS");
                    sor(matrix, b, x, 0.5, 2);
                    system("pause");
                    break;*/
      }
      return 0;
}

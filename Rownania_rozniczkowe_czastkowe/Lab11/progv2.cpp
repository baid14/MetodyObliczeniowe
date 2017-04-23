#include <cstdlib>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
//#include "calerf.h"

using namespace std;

double TMAX = 2.0;
double r = 1.0;
double a = 10.0;
double D = 1.0;

double ppp = r;			//poczatek przedzialu przestrzennego
double kpp = r+a;		//koniec przedzialu przestrzennego
double ppc = 0.0;		//poczatek przedzialu czasowego
double kpc = TMAX;		//koniec przedzialu czasowego

double *przestrzenny;	//krok przestrzenny
double *czasowy;		//krok czasowy

double** CN;
double** KMB;

int N, M;				//ilosci podprzedzialow

double rozwiazanieAnalityczne(int t, int x)
{
	return 1.0 - (r/przestrzenny[x]) * erfc((przestrzenny[x]-r)/(2.0*sqrt(D*czasowy[t])));
}

void thomas1(double *l, double *d, double *u)
{
	for( int i = 1; i < N; i++)
		d[i] -= ( l[i-1] / d[i-1] ) * u[i-1];
}

void thomas2(double *l, double *d, double *u, double *x, double *b)
{
	for( int i = 1; i < N; i++)
		b[i] -= ( l[i-1] / d[i-1] ) * b[i-1];

	x[N-1] = b[N-1] / d[N-1];

	for( int i = N-2; i >= 0; i--)
		x[i] = ( b[i] - u[i] * x[i+1] ) / d[i];
}


void metodaCrankaNicolson(int n)		//N - liczba podprzedzia³ów przedzia³u przestrzennego
{
	N = n;
	double lambda = 1.0;
	double h = 0.1;
	double dt = lambda*h*h;
	M =  ((kpc-ppc)/dt)+1;		//M - liczba podprzedzia³ów przedzia³u czasowego

	przestrzenny = new double[N];	//stworzenie przedzia³u przestrzennego
	for(int i = 0; i < N; i++)
		przestrzenny[i] = r + i*h;

	czasowy = new double[M];		//stworzenie przedzia³u czasowego
	for(int i = 0; i < M; i++)
		czasowy[i] = i*dt;

	CN = new double*[M];
	for (int i = 0; i < M; i++)
		CN[i] = new double[N];

	double *b;
	b = new double[N];					//wektor wyrazow wolnych w ukladzie rownan Ax = b;

	double *l;						   //wektor zawierajacy wartosci podprzekatnej macierzy A
	l = new double[N-1];
	double *d;						   //wektor zawierajacy wartosci przekatnej macierzy A
	d = new double[N];
	double *u;						   //wektor zawierajacy wartosci nadprzekatnej macierzy A
	u = new double[N-1];

	for( int i = 0; i < N; i++)			//warunek poczatkowy
			CN[0][i] = 1.0;

										//dyskretyzacja
	for( int i = 0; i < N-2; i++ )      // wektor l
         l[i] = lambda/2.0 -(lambda/2.0)/(i+r/h);
    l[N-2]=0.0;

    for( int i = 1; i < N-1; i++ )      // wektor d
         d[i] = -(1.0 + lambda);
    d[N-1] = 1.0;
    d[0] = 1.0;

    u[0]=0.0;
	for( int i = 1; i < N-1; i++ )        // wektor u
         u[i] = lambda/2.0 + (lambda/2.0)/(i-1+(r/h));

    thomas1(l, d, u);

    for( int j = 1; j < M; j++ )
	{
		for(int i = 1; i < N-1; i ++)
			b[i] = ((-lambda/2.0) + (lambda/2.0)/(i+r/h))* CN[j-1][i-1] + (lambda - 1.0)*CN[j-1][i] - ((lambda/2.0) + (lambda/2.0)/(i+r/h))*CN[j-1][i+1];
		b[0] = 0.0;															//warunek brzegowy poczatku przedzialu
		b[N-1] = 1.0 - (r/(r+a) * erfc(a/(2.0*sqrt(D*czasowy[j]))));		//warunek brzegowy konca przedzialu

		thomas2(l, d, u, CN[j], b);
	}
	/* pokazywanie macierzy
	for(int i=0; i<M; i++){
		for(int j=0; j<N; j++){
			cout << rozwiazanieAnalityczne(i,j) << " ";
		}
		cout << endl;
	}
	for(int i=0; i<M; i++){
		for(int j=0; j<N; j++){
			cout << CN[i][j] << " ";
		}
		cout << endl;
	}
	*/
	ofstream plik1("anacn.txt");
	for(int i=0; i<N; i++){
		plik1 << przestrzenny[i] << " " << rozwiazanieAnalityczne(M-1,i) << endl;
	}

	ofstream plik("cn.txt");
	for(int i=0; i<N; i++){
		plik << przestrzenny[i] << " " << CN[M-1][i] << endl;
	}
}

void klasycznaMetodaBezposrednia(int n)
{
	N = n;
	double lambda = 0.4;
	double h = 0.1;
	double dt = lambda*h*h;
	M =  ((kpc-ppc)/dt)+1;		//M - liczba podprzedzia³ów przedzia³u czasowego

	przestrzenny = new double[N];	//stworzenie przedzia³u przestrzennego
	for(int i = 0; i < N; i++)
		przestrzenny[i] = r + i*h;

	czasowy = new double[M];		//stworzenie przedzia³u czasowego
	for(int i = 0; i < M; i++)
		czasowy[i] = i*dt;

	KMB = new double*[M];
	for (int i = 0; i < M; i++)
		KMB[i] = new double[N];

	for( int i = 0; i < N; i++)			//warunek poczatkowy
			KMB[0][i] = 1.0;
    double dh;
	for( int j = 1; j < M; j++ )
	{
		KMB[j][0] = 0.0;
		KMB[j][N-1] = 1.0 - (r/(r+a) * erfc(a/(2.0*sqrt(czasowy[j]))));
		dh = r + h;
		for(int i = 1; i < N-1; i++)
        {

            KMB[j][i] = lambda *(1.0 - h / dh)*KMB[j-1][i-1] + (1.0 - ( 2.0 * lambda ))*KMB[j-1][i] + lambda *(( 1.0 + h / dh ))*KMB[j-1][i+1];
        dh +=h;
        }
			}
	/* pokazywanie macierz
	for(int i=0; i<M; i++){
		for(int j=0; j<N; j++){
			cout << rozwiazanieAnalityczne(i,j) << " ";
		}
		cout << endl;
	}
	for(int i=0; i<M; i++){
		for(int j=0; j<N; j++){
			cout << KMB[i][j] << " ";
		}
		cout << endl;
	}
	*/
	ofstream plik1("anakmb.txt");
	for(int i=0; i<N; i++){
		plik1 << przestrzenny[i] << " " << rozwiazanieAnalityczne(M-1,i) << endl;
	}

	ofstream plik("kmb.txt");
	for(int i=0; i<N; i++){
		plik << przestrzenny[i] << " " << KMB[M-1][i] << endl;
	}
}

int main()
{
	metodaCrankaNicolson(100);
	cout << endl;
	klasycznaMetodaBezposrednia(100);
}

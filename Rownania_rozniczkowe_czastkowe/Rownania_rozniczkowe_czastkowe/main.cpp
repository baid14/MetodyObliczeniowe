#include <iostream>
#include <Windows.h>
#include <math.h>
#include "calerf.h"
#include <iomanip> 
#include <fstream>

using namespace std;
//Stałe Wspolczynniki
const double t_min = 0.0, t_max = 2.0;
const double r = 1.0, a = 10.0;
const double D = 1.0, LAMBDA_KMB = 0.4, LAMBDA_LAASONEN = 1.0;
// t_size il poziomow czasowych dla laasonen i kmb 
int t_size_laasonen, t_size_kmb, x_size; 
//pliki do zapisania danych do wykresu
fstream dane_wykres_analityczny, dane_blad_tmax, dane_blad_t_laasonen, dane_blad_kmb;

//Funkcja obliczajaca wartosci analityczne dla Klasycznej metody bezposredeniej
void rozwiazanie_analityczne_K( double **A, double dt , double h)
{
	double tt = dt, dh;
	for ( int i = 1; i < t_size_kmb; i++ )
	{
		dh = 1;
		for ( int j = 0; j < x_size; j++ )
		{
			A[i][j] = 1.0 - ( r / dh ) * calerf_erfc( ( dh - r ) / ( 2.0 * sqrt( D * tt ) ) );
			dh += h;
		}
		tt += dt;
	}
}
//Funkcja obliczajaca wartosci analityczne dla Metody posredniej laasonen
void rozwiazanie_analityczne_L( double **A, double dt, double h )
{
	double tt = dt, dh;
	for ( int i = 1; i < t_size_laasonen; i++ )
	{
		dh = 1;
		for ( int j = 0; j < x_size; j++ )
		{
			A[i][j] = 1.0 - ( r /dh ) * calerf_erfc( ( dh - r ) / ( 2.0 * sqrt( D * tt ) ) );
			dh += h;
		}
		tt += dt;
	}
}

void thomas_redukcja(double *l, double *d, double *u)
{
	for ( int i = 0; i < x_size; i++ )
	{
		d[i + 1] = d[i + 1] - l[i] * ( u[i] / d[i] );
		l[i] = l[i] / d[i];
	}
}

void thomas_rozwiazanie_ukladu( double *l, double *d, double *u, double *b, double *x )
{
	for ( int i = 1; i < x_size ; i++ )
		b[i] = b[i] - l[i - 1] * b[i - 1];

	x[x_size - 1] = b[x_size - 1] / d[x_size - 1];
	for ( int i = x_size - 2; i >= 0; i-- )
		x[i] = ( b[i] - u[i] * x[i + 1] ) / d[i];
}
//obliczenie odpowienich wartosc do macierzy trojdiagonalnej - dyskretyzacja
void dyskretyzacja_Laasonen( double *l, double *d, double *u, double h )
{
	double dh = r + h;

	d[0] = 1.0;
	u[0] = 0.0;
	for ( int i = 1; i <= x_size - 2; i++ )
	{
		l[i - 1] = LAMBDA_LAASONEN * ( 1.0 - ( h / dh ) ); //wektor l
		d[i] = -( 1.0 + ( 2 * LAMBDA_LAASONEN ) ); //wektor d
		u[i] = LAMBDA_LAASONEN * ( 1.0 + h / dh ); //wektor u
		dh += h;
	}
	d[x_size - 1] = 1.0;
	l[x_size - 2] = 0.0;		
}
// Funkcja realizujaca Metoda bezposrednia Laasonen
void metoda_posrednia_laasonen( double ** A, double h )
{
	double *l, *d, *u, *b, *x;
	l = new double[x_size - 1];
	d = new double[x_size];
	u = new double[x_size - 1];
	b = new double[x_size];
	x = new double[x_size];

	dyskretyzacja_Laasonen(l,d,u, h);
	thomas_redukcja( l, d, u );
	for ( int i = 1; i < t_size_laasonen; i++ )
	{
		for ( int j = 1; j < x_size - 1; j++ )
			b[j] = -A[i - 1][j]; // wektor b

		b[0] = A[i][0];
		b[x_size - 1] = A[i][x_size - 1];
		thomas_rozwiazanie_ukladu( l, d, u, b, x );
		for ( int k = 1; k < x_size - 1 ; k++ )
			A[i][k] = x[k];	//wypelnianie macierzy aktualnym poziomem czasowym
	}

	delete[] l;
	delete[] d;
	delete[] u;
	delete[] x;
	delete[] b;
}
// Funkcja realizujaca Klasyczna metode bezposrednia
void metoda_KMB( double **A, double h )
{
	double dh;
	for ( int i = 1; i < t_size_kmb; i++ )
	{
		dh = r + h;
		for ( int j = 1; j < x_size - 1; j++ )
		{ //obliczanie nowego poziomu czasowego zgodnie z wzorem
			A[i][j] = A[i - 1][j - 1] * LAMBDA_KMB * ( 1.0 - h / dh ) 
				+ A[i - 1][j] * ( 1.0 - ( 2.0 * LAMBDA_KMB ) ) 
				+ A[i - 1][j + 1] * LAMBDA_KMB * ( 1.0 + h / dh );
			dh += h;
		}
	}
}
//warynek brzegowy
double prawy( double t )
{
	return 1.0 - ( r / ( r + a ) ) * calerf_erfc( a / ( 2.0 * sqrt( D * t ) ) );
}
//Funkcja realizujaca zadanie
void rozwiazanie()
{
	dane_blad_tmax.open( "Dane_blad_tmax.txt", fstream::out );

	double blad_laasonen, blad_kmb, max_blad_lassonen, max_blad_kmb, h = 0.5;
 
	for ( h; h > 0.01; h -= 0.01 )
	{
		double dt_laasonen = ( LAMBDA_LAASONEN * h * h ) / D;

		x_size = a / h;

		//Metoda_lasonen_dane
		double **Analityczne_L, **Laasonen;

		t_size_laasonen = t_max / dt_laasonen ;

		Analityczne_L = new double *[t_size_laasonen];
		Laasonen = new double *[t_size_laasonen];

		for ( int i = 0; i < t_size_laasonen; i++ )
		{
			Analityczne_L[i] = new double[x_size];
			Laasonen[i] = new double[x_size];
		}

		//Uzupelnianie tablic warunkami poczatkowymi
		//Warunek poczatkowy
		for ( int i = 0; i < x_size; i++ )
		{
			Laasonen[0][i] = 1.0;
		}
		//1 Warunek brzegowy 
		for ( int i = 1; i < t_size_laasonen; i++ )
		{
			Laasonen[i][0] = 0.0;
		}
		//2 Warunek brzegowy
		double tt = dt_laasonen;
		for ( int i = 1; i < t_size_laasonen; i++ )
		{
			Laasonen[i][x_size - 1] = prawy( tt );
			tt += dt_laasonen;
		}

		//Metoda_KMB
		double **Analityczne_K, **KMB;
		double dt_kmb = ( LAMBDA_KMB * h *h ) / D;

		t_size_kmb = t_max / dt_kmb ;

		Analityczne_K = new double *[t_size_kmb];
		KMB = new double *[t_size_kmb];

		for ( int i = 0; i < t_size_kmb; i++ )
		{
			Analityczne_K[i] = new double[x_size];
			KMB[i] = new double[x_size];
		}

		//Uzupelnianie tablic warunkami poczatkowymi
		//Warunek poczatkowy
		for ( int i = 0; i < x_size; i++ )
		{
			KMB[0][i] = 1.0;
		}
		//1 Warunek brzegowy 
		for ( int i = 1; i < t_size_kmb; i++ )
		{
			KMB[i][0] = 0.0;
		}
		//2 Warunek brzegowy
		tt = dt_kmb;
		for ( int i = 1; i < t_size_kmb; i++ )
		{
			KMB[i][x_size - 1] = prawy( tt );
			tt += dt_kmb;
		}
		//wpelnienie tablic poziomami czasowymi
		rozwiazanie_analityczne_L( Analityczne_L, dt_laasonen, h);
		metoda_posrednia_laasonen( Laasonen, h );
		rozwiazanie_analityczne_K( Analityczne_K, dt_kmb, h );
		metoda_KMB( KMB, h );
		cout << h << endl;
		//gdy h == 0.1 zapisuje obliczenia do wykresu funkcji 
		// oraz licze maksymalna wartosc bezwgledna bledu w funkcji czasu t
		if ( h < 0.101 && h > 0.099 )
		{
			cout << endl << t_size_laasonen << " " << t_size_kmb << endl;
			cout << endl << "Rozwiazania" << endl;
			cout << "  x | Analityczne_L |      Laasonen | Analityczne_K |         KMB | " << endl;
			cout << "-------------------------------------------------------------------" << endl;
			int index = 150;
			double dh = r + h;
			dane_wykres_analityczny.open( "dane_wykres_analityczny.txt", fstream::out );
			for ( int i = 0; i < x_size; i++ )
			{
				dane_wykres_analityczny << dh << " " << Analityczne_L[index][i] << " " << Laasonen[index][i] << " " <<
					Analityczne_K[index][i] << " " << KMB[index][i] << endl;
				dh += h;
				cout.width( 4 );
				cout << dh << "|";
				cout.width( 15 );
				cout << setprecision( 10 ) << Analityczne_L[index][i] << "|";
				cout.width( 15 );
				cout << setprecision( 10 ) << Laasonen[index][i] << "|";
				cout.width( 15 );
				cout << setprecision( 10 ) << Analityczne_K[index][i] << "|";
				cout.width( 15 );
				cout << setprecision( 10 ) << KMB[index][i] << "|" << endl;
			}
			dane_wykres_analityczny.close();

			double blad_laasonen, blad_kmb, max_blad_lassonen, max_blad_kmb;
			
			dane_blad_t_laasonen.open( "dane_blad_t_laasonen.txt", fstream::out );
			for ( int i = 1; i < t_size_laasonen; i++ )
			{
				max_blad_lassonen = 0.0;
				for ( int j = 0; j < x_size; j++ )
				{
					blad_laasonen = fabs( Analityczne_L[i ][j] - Laasonen[i][j] );
					if ( blad_laasonen > max_blad_lassonen )
						max_blad_lassonen = blad_laasonen;
				}
				dane_blad_t_laasonen << i * dt_laasonen << " " << max_blad_lassonen << endl;
			}
			dane_blad_t_laasonen.close();

			dane_blad_kmb.open( "dane_blad_kmb.txt", fstream::out );
			for ( int i = 1; i < t_size_kmb; i++ )
			{
				max_blad_kmb = 0.0;
				for ( int j = 0; j < x_size; j++ )
				{
					blad_kmb = fabs( Analityczne_K[i][j] - KMB[i][j] );
					if ( blad_kmb > max_blad_kmb )
						max_blad_kmb = blad_kmb;
				}
				dane_blad_kmb << i*dt_kmb << " " << max_blad_kmb << endl;
			}
		}
		max_blad_lassonen = 0.0;
		max_blad_kmb = 0.0;
		//Maksymalna wartosc bezwzglednego bledu obserwacji dla t_max w funkcji kroku przestrzennnego h
		for ( int i = 0; i < x_size; i++ )
		{
			blad_laasonen = fabs( Analityczne_L[t_size_laasonen - 1][i] - Laasonen[t_size_laasonen - 1][i] );
			if ( blad_laasonen > max_blad_lassonen )
				max_blad_lassonen = blad_laasonen;

			blad_kmb = fabs( Analityczne_K[t_size_kmb - 1][i] - KMB[t_size_kmb - 1][i] );
			if ( blad_kmb > max_blad_kmb )
				max_blad_kmb = blad_kmb;

		}
		dane_blad_tmax << log10( h ) << " " << log10( max_blad_lassonen ) << " " << log10( max_blad_kmb ) << endl;

		//usuwanie tablic z pamieci
		for ( int i = 0; i < t_size_kmb; i++ )
		{
			delete[] Analityczne_K[i];
			delete[] KMB[i];
		}

		for ( int i = 0; i < t_size_laasonen; i++ )
		{
			delete[] Analityczne_L[i];
			delete[] Laasonen[i];
		}

		delete[] Analityczne_K;
		delete[] Analityczne_L;
		delete[] Laasonen;
		delete[] KMB;		

	}
;}

int main()
{
	rozwiazanie();
	system( "Pause" );
	return 0;
}
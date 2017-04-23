#define _USE_MATH_DEFINES

#include <iostream>
#include <Windows.h>
#include <fstream>
#include <math.h>

using namespace std;

const int il_wezlow = 15;
const double X_MIN = -1.0, X_MAX = 1.0;
double const skok_il_wezlow = ( X_MAX - X_MIN ) / ( il_wezlow - 1 );

double Funkcja( double x)
{
	return x / ( 1.0 + 25.0 * pow( fabs( x ), 3 ) );
}

double Lagrange( double x, double *xi, double *yi )
{
	double iloczyn[il_wezlow];
	double wynik = 0.0;

	for ( int i = 0 ; i < il_wezlow ; i++ )
	{
		iloczyn[i] = 1.0;
		for ( int j = 0 ; j < il_wezlow ; j++ )
		if ( j != i )
			iloczyn[i] = iloczyn[i] * ( ( x - xi[j] ) / ( xi[i] - xi[j] ) );
	}

	for ( int i = 0 ; i < il_wezlow ; i++ )
		wynik = wynik + yi[i] * iloczyn[i];

	return wynik;
}

double Czybyszew( double x )
{
	double xCz[il_wezlow];
	double yCz[il_wezlow];

	for ( int i = 0; i < il_wezlow; i++ )
	{
		xCz[i] = cos( ( ( 2.0 * i + 1.0 ) * M_PI ) / ( 2.0 * il_wezlow ) ); // wyznaczanie kolejnego wezla
		yCz[i] = Funkcja( xCz[i] );
	}

	return Lagrange( x, xCz, yCz );

}

double Rownoodlegle( double x )
{
	double xi[il_wezlow];
	double yi[il_wezlow];

	int i = 0;

	for ( double z = X_MIN; z <= X_MAX; z += skok_il_wezlow )
	{
		xi[i] = z; // poczatek przedzialu
		yi[i] = Funkcja( z ); // wyliczenie funkcji dla wezla 
		i++;
	}
	return Lagrange( x, xi, yi );
}

void rozwiazanie()
{
	fstream dane;

	const int N = 100; 

	double x = X_MIN;
	double w_Funkcja[N + 1], w_Czybyszew[N+1], w_Rownoodlegle[N+1];
	double skok = ( X_MAX - X_MIN ) / N;

	cout << "  x | Analityczne |      Czybyszew | Rownoodlegle | " << endl;
	cout << "----------------------------------------------------" << endl;

	dane.open( "wyniki.txt", fstream::out );
	for ( int i = 0 ; i < N + 1; i++ )
	{
		w_Funkcja[i] = Funkcja( x );
		w_Czybyszew[i] = Czybyszew(x);
		w_Rownoodlegle[i] = Rownoodlegle( x );
		dane << x << " " << w_Funkcja[i] << " " << w_Czybyszew[i] << " " << w_Rownoodlegle[i] << endl;
		cout.width( 4 );
		cout << x << "|";
		cout.width( 15 );
		cout << w_Funkcja[i] << "|";
		cout.width( 15 );
		cout << w_Czybyszew[i] << "|";
		cout.width( 15 );
		cout << w_Rownoodlegle[i] << "|" << endl;
		x += skok;
	}
	dane.close();

}

int main()
{
	rozwiazanie();
	system( "Pause" );
	return 0;
}
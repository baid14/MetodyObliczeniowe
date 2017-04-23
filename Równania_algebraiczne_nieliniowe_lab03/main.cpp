#include <iostream>
#include <Windows.h>
#include <math.h>

using namespace std;

const int il_petli = 100; 
const double eps = 1e-7;

double pierwsza_funkcja( double x )
{
	return pow( cos( x / 4 ), 2 ) - x;
}
double pierwsza_funkcja_pochodna( double x )
{
	return -0.25 * sin( x / 4.0 ) - 1;
}
double pierwsza_funckja_picard( double x )
{
	return pow( cos( x / 4 ), 2 );
}
double pierwsza_funkcja_picard_pochodna( double x )
{
	return -0.25 * sin( x / 4.0 );
}


double druga_funkcja( double x )
{
	return exp( x ) - exp( -x ) + x - 1;
}
double druga_funkcja_pochodna( double x )
{
	return exp( x ) + exp( x ) + 1;
}
double druga_funkcja_picard( double x )
{
	return exp(-x) - exp(x) + 1;
}
double druga_funkcja_picard_pochodna( double x)
{
	return -exp( -x ) - exp( x );
}

//Metoda iteracji prostych - Picarda
void metoda_Picarda( double (*funkcja) (double), double( *funkcja_picard )( double ), double( *funkcja_pochodna )( double ), double x0 )
{
	cout << endl << "Metoda Picarda:" << endl;

	if ( fabs( funkcja_pochodna( x0 ) ) >= 1 )
	{
		cout << "Funkcja jest rozbiezna." << endl;
	} else
	{
		double x1, fx0;
		double EST;
		cout << "  n |            xn | EST(xn - x*)  |       |f(xn)| |" << endl;
		cout << "-----------------------------------------------------" << endl;
		for ( int i = 0; i < il_petli; i++ )
		{
			fx0 = fabs(funkcja( x0 ));
			x1 = funkcja_picard( x0 );
			EST = fabs( x1 - x0 );
			x0 = x1;
			
			cout.width( 4 );
			cout << i << "|";
			cout.width( 15 );
			cout << x0 << "|";
			cout.width( 15 );
			cout << EST << "|";
			cout.width( 15 );
			cout << fx0 << "|" << endl;

			if ( EST < eps ) break;
			if ( fx0 < eps ) break;
		}
	}

}

//Metoda Bisekcji
void metoda_Bisekcji(double (*funkcja)(double), double a, double b)
{
	cout << endl << "Metoda Bisekcji:" << endl;

	if ( funkcja( a ) * funkcja( b ) > 0 )
	{
		cout << "Podano zly przedzial." << endl;
	} else
	{
		double EST, xn, fxn;

		cout << "  n |            xn | EST(xn - x*)  |       |f(xn)| |" << endl;
		cout << "-----------------------------------------------------" << endl;

		for ( int i = 0; i < il_petli; i++ )
		{
			xn = ( a + b ) / 2;
			EST = fabs( b - a ) / 2 ;
			fxn = funkcja( xn );
			if ( funkcja( a ) * fxn < 0 )
				b = xn;
			else
				a = xn;

			cout.width( 4 );
			cout << i << "|";
			cout.width( 15 );
			cout << xn << "|";
			cout.width( 15 );
			cout << EST << "|";
			cout.width( 15 );
			cout << fabs(fxn) << "|" << endl;

			if ( EST < eps ) break;
			if ( fabs(fxn) < eps ) break;
		}
	}

}

//Metoda Newtona
void metoda_Newtona( double( *funkcja )( double ), double(*funkcja_pochodna)(double), double xn )
{
	cout << endl << "Metoda Newtona:" << endl;

	double fxn, xn1, EST;

	cout << "  n |            xn | EST(xn - x*)  |       |f(xn)| |" << endl;
	cout << "-----------------------------------------------------" << endl;

	for ( int i = 0; i < il_petli; i++ )
	{
		xn1 = xn;
		fxn = funkcja( xn1 );
		xn = xn1 - fxn / funkcja_pochodna( xn1 );

		EST = fabs( xn1 - xn );

		cout.width( 4 );
		cout << i << "|";
		cout.width( 15 );
		cout << xn << "|";
		cout.width( 15 );
		cout << EST << "|";
		cout.width( 15 );
		cout << fabs(fxn) << "|" << endl;

		if ( EST < eps ) break;
		if ( fabs( fxn ) < eps ) break;
	}
}

//Metoda Siecznych
void metoda_Siecznych( double( *funkcja )( double ), double xn0, double xn1)	//x0 ten dalszy xn1 blizszy x2 obliczany
{
	cout << endl << "Metoda Siecznych:" << endl;

	double xn2, EST, fxn2;

	cout << "  n |            xn | EST(xn - x*)  |       |f(xn)| |" << endl;
	cout << "-----------------------------------------------------" << endl;

	for ( int i = 0; i < il_petli; i++ )
	{
		xn2 = xn0 - funkcja( xn0 ) / ( ( funkcja( xn0 ) - funkcja( xn1 ) ) / ( xn0 - xn1 ) );
		xn1 = xn0; 
		xn0 = xn2;
		fxn2 = fabs(funkcja( xn0 ));
		EST = fabs( xn0 - xn1 );

		cout.width( 4 );
		cout << i << "|";
		cout.width( 15 );
		cout << xn2 << "|";
		cout.width( 15 );
		cout << EST << "|";
		cout.width( 15 );
		cout << fxn2 << "|" << endl;

		if (  EST < eps ) break;
		if ( fxn2  < eps ) break;
	}

}

int main()
{
	cout << "Pierwsza funkcja: cos^2(x/4) - x = 0" << endl;
	metoda_Picarda( pierwsza_funkcja, pierwsza_funckja_picard, pierwsza_funkcja_picard_pochodna, 5 );
	metoda_Bisekcji( pierwsza_funkcja, -5, 5 );
	metoda_Newtona( pierwsza_funkcja, pierwsza_funkcja_pochodna, 5 );
	metoda_Siecznych( pierwsza_funkcja, 5, 10 );
	system( "pause" );

	cout << "Druga funkcja: exp(x) - exp(-x) + x - 1 = 0 " << endl;
	metoda_Picarda( druga_funkcja, druga_funkcja_picard, druga_funkcja_picard_pochodna, 5 );
	metoda_Bisekcji( druga_funkcja, -5, 5 );
	metoda_Newtona( druga_funkcja, druga_funkcja_pochodna, 5 );
	metoda_Siecznych( druga_funkcja, 5, 10 );

	system( "pause" );
	return 0;
}
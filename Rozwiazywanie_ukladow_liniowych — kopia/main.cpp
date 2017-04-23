#include <iostream>
#include <Windows.h>
#include <math.h>

using namespace std;

const int n = 3;
const int il_petli = 100;
const double eps = 1e-7;

void wypiszMac( double **mac )
{

	for ( int i = 0; i<n; i++ )
	{
		cout << "| ";
		for ( int j = 0; j<n; j++ )
		{
			//printf( "%lf, ", mac[indeksy[i]][j] );
			cout.width( 4 );
			cout << mac[i][j];
		}
		cout << " |" << endl;
	}
}

void wypiszWek( double *wek)
{
	for ( int i = 0; i<n; i++ )
	{

		cout << "| ";
		cout.width( 4 );
		//printf( "%lf, ", wek[indeksy[i]] );
		cout << wek[i];
		cout << "| ";
		cout << endl;
	}
}

void wypelnienie( double **A, double *b, double *x )
{
	A[0][0] = 10.0; A[0][1] = 1.0; A[0][2] = 2.0; 
	A[1][0] = 1.0; A[1][1] = 20.0; A[1][2] = 5.0;
	A[2][0] = 3.0; A[2][1] = 4.0; A[2][2] = 30.0; 

	b[0] = 8.0; b[1] = -4.0; b[2] = -27.0;

	x[0] = 2.0; x[1] = 2.0; x[2] = 2.0;
}

double est( double *x, double *x_nowe )
{
	x[0] = fabs( x[0] - x_nowe[0]);
	x[1] = fabs( x[1] - x_nowe[1]);
	x[2] = fabs( x[2] - x_nowe[2]);
	
	double max = x[0];
	if ( x[1] > max ) max = x[1];
	if ( x[2] > max ) max = x[2];

	return max;
}

double residuum( double **A, double *b, double *x_nowe )
{
	double Ax[3];
	
	Ax[0] = fabs( (A[0][0] * x_nowe[0] + A[0][1] * x_nowe[1] + A[0][2] * x_nowe[2]) - b[0] );
	Ax[1] = fabs( (A[1][0] * x_nowe[0] + A[1][1] * x_nowe[1] + A[1][2] * x_nowe[2]) - b[1] );
	Ax[2] = fabs( (A[2][0] * x_nowe[0] + A[2][1] * x_nowe[1] + A[2][2] * x_nowe[2]) - b[2] );

	double max = Ax[0];
	if ( Ax[1] > max ) max = Ax[1];
	if ( Ax[2] > max ) max = Ax[2];
	
	return max;
}

void metoda_Jacobiego( double **A, double *b, double *x )
{
	double *x_nowe = new double[n]; //nowe przyblizenia
	double suma = 0.0;
	double EST = 0.0, RESIDUUM = 0.0;

	cout << endl << endl << " \t Metoda Jacobiego" << endl;
	cout << "  n |            x1 |            x2 |            x3 |            EST |     RESIDIUM |" << endl;
	cout << "----------------------------------------------------------------------------------" << endl;

	for ( int iter = 0; iter < il_petli; iter++ )
	{
		for ( int i = 0; i < n; i++ )
		{
			suma = 0.0;
			for ( int j = 0; j < n; j++ )
				if ( j != i )
					suma += A[i][j] * x[j];
						
			x_nowe[i] = ( 1.0 / A[i][i] ) * ( b[i] - suma );		
		}

		EST = est( x, x_nowe );
		RESIDUUM = residuum( A, b, x_nowe );

		for ( int i = 0; i < n; i++ )
			x[i] = x_nowe[i];

		cout.width( 4 );
		cout << iter << "|";
		cout.width( 15 );
		cout << x_nowe[0] << "|";
		cout.width( 15 );
		cout << x_nowe[1] << "|";
		cout.width( 15 );
		cout << x_nowe[2] << "|";
		cout.width( 15 );
		cout << EST << "|";
		cout.width( 15 );
		cout << RESIDUUM << "|" << endl;

		if ( EST < eps )	
			break;
			
		if ( RESIDUUM < eps )
			break;
		}
	cout << "----------------------------------------------------------------------------------" << endl;
	}

void metoda_Gaussa_Seidela( double **A, double *b, double *x )
{
	double *x_poprz = new double[n]; //stare wart
	double suma = 0.0;
	double EST = 0.0, RESIDUUM = 0.0;

	cout << endl << endl << "\t Metoda Gaussa_Seidela" << endl;
	cout << "  n |            x1 |            x2 |            x3 |            EST |     RESIDIUM |" << endl;
	cout << "----------------------------------------------------------------------------------" << endl;

	for ( int iter = 0; iter < il_petli; iter++ )
	{
		for ( int i = 0; i < n; i++ )
		{
			suma = 0.0;
			for ( int j = 0; j < n; j++ )
				if ( j != i )
					suma += A[i][j] * x[j];

			x_poprz[i] = x[i];
			x[i] = ( 1.0 / A[i][i] ) * ( b[i] - suma );
		}

		EST = est( x_poprz, x );
		RESIDUUM = residuum( A, b, x );


		cout.width( 4 );
		cout << iter << "|";
		cout.width( 15 );
		cout << x[0] << "|";
		cout.width( 15 );
		cout << x[1] << "|";
		cout.width( 15 );
		cout << x[2] << "|";
		cout.width( 15 );
		cout << EST << "|";
		cout.width( 15 );
		cout << RESIDUUM << "|" << endl;

		if ( EST < eps )
			break;

		if ( RESIDUUM < eps )
			break;
	}
	cout << "----------------------------------------------------------------------------------" << endl;
}

void metoda_SOR( double **A, double *b, double *x )
{
	double *x_nowe = new double[n]; //nowe przyblizenia
	double *x_poprz = new double[n]; //stare wartosc
	double suma = 0.0, omega = 0.5;
	double EST = 0.0, RESIDUUM = 0.0;

	cout << endl << endl << "\t Metoda SOR" << endl;

	cout << "  n |            x1 |            x2 |            x3 |            EST |     RESIDIUM |" << endl;
	cout << "----------------------------------------------------------------------------------" << endl;

	for ( int iter = 0; iter < il_petli; iter++ )
	{
		for ( int i = 0; i < n; i++ )
		{
			suma = 0.0;
			for ( int j = 0; j < n; j++ )
			if ( j != i )
				suma += A[i][j] * x[j];

			x_poprz[i] = x[i];
			x_nowe[i] = ( 1.0 - omega ) * x[i] + ( omega / A[i][i] ) * ( b[i] - suma );
			x[i] = x_nowe[i];
		}

		EST = est( x_poprz, x_nowe );
		RESIDUUM = residuum( A, b, x_nowe );


		cout.width( 4 );
		cout << iter << "|";
		cout.width( 15 );
		cout << x[0] << "|";
		cout.width( 15 );
		cout << x[1] << "|";
		cout.width( 15 );
		cout << x[2] << "|";
		cout.width( 15 );
		cout << EST << "|";
		cout.width( 15 );
		cout << RESIDUUM << "|" << endl;

		if ( EST < eps )
			break;

		if ( RESIDUUM < eps )
			break;
	}
	cout << "----------------------------------------------------------------------------------" << endl;
}

void rozwiazanie()
{
	double **A, *b, *x;

	A = new double*[n];

	for ( int i = 0; i < n; i++ )
		A[i] = new double[n];

	b = new double[n];

	x = new double[n];

	wypelnienie( A , b , x );

	cout << "Macierz A:" << endl;
	wypiszMac( A );
	cout << endl;
	cout << " Wektor b:" << endl;
	wypiszWek( b );
	cout << " Wektor x:" << endl;
	wypiszWek( x );

	metoda_Jacobiego( A, b, x );
	wypelnienie( A, b, x );
	metoda_Gaussa_Seidela( A, b, x );
	wypelnienie( A, b, x );
	metoda_SOR( A, b, x );

	for ( int i = 0; i < n; i++ )
	{
		delete[] A[i];
	}

	delete[] A;
	delete[] b;
	delete[] x;

}

int main()
{
	rozwiazanie();
	system( "pause" );
	return 0;
}
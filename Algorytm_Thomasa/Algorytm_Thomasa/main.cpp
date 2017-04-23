#include <iostream>
#include <Windows.h>

using namespace std;

void wypelnienie( double *l, double *d, double *u, double *b )
{
	l[0] = 2.0 / 3.0; l[1] = 4.0 / 5.0; l[2] = 6.0 / 7.0; l[3] = 8.0 / 9.0; l[4] = 10.0 / 11.0;

	d[0] = 1.0; d[1] = 2.0; d[2] = 3.0; d[3] = 4.0; d[4] = 5.0; d[5] = 6.0;

	u[0] = 1.0 / 2.0; u[1] = 3.0 / 4.0; u[2] = 5.0 / 6.0; u[3] = 7.0 / 8.0; u[4] = 9.0 / 10.0;

	b[0] = 7.0 / 2.0; b[1] = -3.0; b[2] = 11.0 / 2.0; b[3] = -191.0 / 28.0; b[4] = 193.0 / 30.0; b[5] = -46.0 / 11.0;
}

void wypiszWek( double *wek, int n)
{
	for ( int i = 0; i < n; i++ )
	{
		cout << "| ";
		cout.width( 4 );
		//printf( "%lf, ", wek[i] );
		cout << wek[i];
		cout.width( 4 );
		cout << "| ";
		cout << endl;
	}
}

void thomas( double *l, double *d, double *u, int n )
{
	for ( int i = 0; i < n; i++ )
	{
		d[i+ 1] = d[i + 1] - l[i] * ( u[i] / d[i] );
		l[i] = l[i] / d[i];
	}
}

void rozwiazanie_ukladu( double *l, double *d, double *u, double *b, double *x, int n )
{
	for ( int i = 1; i < n + 1 ; i++ )
		b[i] = b[i] - l[i - 1] * b[i - 1];

	x[n] = b[n ] / d[n ];
	for ( int i = n -1; i >= 0; i-- )
	{

		x[i] = ( b[i] - u[i] * x[i + 1] ) / d[i];
	}
}

void roziwazanie()
{
	double *l, *d, *u, *b, *x;
	int n = 5;

	l = new double[n];
	d = new double[n + 1];
	u = new double[n];
	b = new double[n + 1];
	x = new double[n + 1];
	
	wypelnienie( l, d, u, b );
	
	cout << endl << "Wektor l " << endl;
	wypiszWek( l, n);
	cout << endl << "Wektor d " << endl;
	wypiszWek( d, n + 1 );
	cout << endl << "Wektor u " << endl;
	wypiszWek( u, n );
	cout << endl << "Wektor b " << endl;
	wypiszWek( b, n + 1 );
	cout << endl;

	thomas( l, d, u, n );
	cout << endl << "Wektor l " << endl;
	wypiszWek( l, n );
	cout << endl << "Wektor d " << endl;
	wypiszWek( d, n + 1 );
	cout << endl << "Wektor u " << endl;
	wypiszWek( u, n );

	rozwiazanie_ukladu( l, d, u, b, x, n );

	cout << endl << "Wektor x " << endl;
	wypiszWek( x, n + 1 );
	cout << endl;

	delete[] l;
	delete[] d;
	delete[] u;
	delete[] x;
	delete[] b;

}


int main()
{
	roziwazanie( );
	system( "pause" );
	return 0;
}
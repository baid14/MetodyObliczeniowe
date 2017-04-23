#include <iostream>
#include <Windows.h>
#include <math.h>

using namespace std;

const double eps = 1e-7;

void wypiszMac( double **mac, int indeksy[4] )
{
	for ( int i = 0; i<4; i++ )
	{
		cout << "| ";
		for ( int j = 0; j<4; j++ )
		{
			cout.width( 4 );
			cout << mac[indeksy[i]][j];
		}
		cout << " |" << endl;
	}
}

void wypiszWek( double *wek, int indeksy[4] )
{
	for ( int i = 0; i<4; i++ )
	{
		cout << "| ";
		cout.width( 4 );
		cout << wek[indeksy[i]] << " |" << endl;
	}
}

void wypelnienie( double **A , double *b)
{
	A[0][0] = 2.0; A[0][1] = 6.0; A[0][2] = 8.0; A[0][3] = 4.0;
	A[1][0] = 1.0; A[1][1] = 3.0; A[1][2] = 5.0; A[1][3] = 7.0;
	A[2][0] = 4.0; A[2][1] = 2.0; A[2][2] = 6.0; A[2][3] = 8.0;
	A[3][0] = 3.0; A[3][1] = 5.0; A[3][2] = 7.0; A[3][3] = 1.0;

	b[0] = -2; b[1] = -18; b[2] = -14; b[3] = 10;
}

bool czesciowy_wybor_elementu_podstawowego( double **A, int indeks[4], int x )
{
	double max = 0;
	int indeks_max = x, temp;
	bool znaleziono = false;

	for ( int i = x; i < 4; ++i )
	{
		if ( fabs( A[indeks[i]][x] ) > max )
		{
			max = fabs( A[indeks[i]][x] );
			indeks_max = i;
			znaleziono = true;
		}
	}

	if ( znaleziono )
	{
		temp = indeks[indeks_max];
		indeks[indeks_max] = indeks[x];
		indeks[x] = temp;
	}

	return znaleziono;

}
void dekompozycjaLU()
{
	double **A, **L;
	double *b, y[4], x[4], wspl, temp;
	int indeksy[4];
	int n = 4, i, j,k;


	A = new double *[n];
	L = new double *[n];
	b = new double[n];
	//x = new double[n];
	//x[0] = 1; x[1] = 1; x[2] = 1; x[3] =1;

	for ( i = 0; i < n; ++i )
		indeksy[i] = i;

	for ( i = 0; i < n; ++i )
		A[i] = new double[n];

	for ( i = 0; i < n; ++i )
		L[i] = new double[n];

	for ( i = 0; i < n; ++i )
	{
		for ( j = 0; j < n; ++j )
			L[i][j] = 0;
	}

	wypelnienie( A, b );
	cout << endl;
	wypiszMac( A, indeksy );
	cout << endl;

	for (  i = 0; i < 3; i++ )
	{
		for (  j = i + 1; j < 4; j++ )
		{
			if ( A[indeksy[i]][indeksy[i]] == 0 )
			{
				if ( !czesciowy_wybor_elementu_podstawowego( A, indeksy, i )  )
				{
					break;
				}
			}

			L[j][i] = A[indeksy[j]][i] / A[indeksy[i]][i];
			//wypiszMac( L, indeksy );
			//cout << L[j][i] << endl;
			wspl = A[indeksy[j]][i] / A[indeksy[i]][i];
			for (  k = i; k < 4; k++ )
			{
				A[indeksy[j]][k] = A[indeksy[j]][k] - A[indeksy[i]][k] * wspl; //U
			}
		}
	}

	for ( i = 0; i < 4; ++i )
	{
		L[indeksy[i]][i] = 1;
	}

	cout << " Macierz U" << endl;
	wypiszMac( A, indeksy );
	cout << endl;
	cout << " Macierz L" << endl;
	wypiszMac( L, indeksy );
	cout << endl;
	cout << " wektor b" << endl;
	wypiszWek( b, indeksy );
	cout << endl;

	for ( i = 0; i < 4; i++ )
	{
		temp = 0;
		for ( j = 0; j < i; j++ )
		{
			temp += L[indeksy[i]][j] * y[indeksy[j]];
		}
		y[indeksy[i]] = ( b[indeksy[i]] - temp );

	}
	cout << " wektor y" << endl;

	wypiszWek( y, indeksy );
	cout << endl;
	for ( i = 4 - 1; i >= 0; i-- )
	{
		temp = 0;
		for ( j = i + 1; j < 4; j++ )
		{
			temp += A[indeksy[i]][j] * x[indeksy[j]];
			cout << x[indeksy[j]] << endl;
		}
		x[indeksy[i]] = ( y[indeksy[i]] - temp ) / A[indeksy[i]][i];
	}



	cout << endl;
	wypiszWek( x, indeksy );

}

int main()
{
	dekompozycjaLU();
	system( "pause" );
	return 0;
}

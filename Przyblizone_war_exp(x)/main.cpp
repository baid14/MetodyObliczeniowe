#include <iostream>
#include <Windows.h>
#include <math.h>

using namespace std;

double szereg_taylora( int x, int n )
{
	double suma = 1.0;
	double silnia = 1.0;
	double potega = 1.0;
	
	for ( int i = 1; i <= n; ++i )
	{
		silnia *= i;
		potega *= x;
		suma += potega / silnia;
	}
	return suma;
	
}

double szereg_taylora_poprawka( double x, int n )
{
	double suma = 1.0;
//	double wyrazy = 1.0;
	double silnia = 1.0;
	double potega = 1.0;

	if ( x >= 0 )
	{
		for ( int i = 1; i <= n; ++i )
		{
			silnia *= i;
			potega *= x;
			suma += potega / silnia;

			//wyrazy *= x / i;
			//suma += wyrazy;
		}
		return suma;
	} else 
		return 1 / szereg_taylora_poprawka( -x, n );
}


int main()
{
	double wd, wp, rerr, wp2, rerr2;
	int n = 180;
	cout << "x   | Wd            | Wp            | Rerr          | Wp~           | Rerr~         " << endl;
	cout << "-----------------------------------------------------------------------------------" << endl;

	for ( int x = -30; x <= 30; ++x )
	{
		wd = exp( x );

		wp = szereg_taylora( x, n );
		rerr = abs( ( wp - wd ) / wd );

		wp2 = szereg_taylora_poprawka( x, n );
		rerr2 = abs( ( wp2 - wd ) / wd );

		cout.width( 4 );
		cout << x << "|";
		cout.width( 15 );
		cout << wd << "|";
		cout.width( 15 );
		cout << wp << "|";
		cout.width( 15 );
		cout << rerr << "|" ;

		cout.width( 15 );
		cout << wp2 << "|";
		cout.width( 15 );
		cout << rerr2 << endl;
	}

	getchar();
	return 0;
}
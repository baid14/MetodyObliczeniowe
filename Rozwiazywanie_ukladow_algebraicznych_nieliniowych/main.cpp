#include <iostream>
#include <math.h>

using namespace std;

const int il_petli = 100;
const double eps = 1e-7;

double funkcja_1( double x, double y, double z )
{
	return x*x + y*y - 1;
}
double funkcja_2( double x, double y, double z )
{
	return 2 * x*x - y + 1/2;
}
double funkcja_3( double x, double y, double z )
{
	return tan( x / 4 ) + y*y*y - z*z*z*z*z;
}

/*
double funkcja_1x( double x )
{
	
}
double funkcja_2x( double x )
{
	return 4 * x;
}
double funkcja_3x( double x )
{
	return 1 + tan( x / 4 )*tan( x / 4 );
}

double funkcja_1y( double y )
{
	return 2y;
}
double funkcja_2y( double y )
{
	return -1;
}
double funkcja_3y( double y )
{
	return 3 * y*y;
}

double funkcja_1z( double z )
{
	return 0;
}
double funkcja_2z( double z )
{
	return 0;
}
double funkcja_3y( double z )
{
	return -5 * z*z*z*z*z;
}
*/

double norma( double x, double y, double z )
{
	double norma = fabs( x );
	if ( fabs( y ) > norma ) norma = fabs( y );
	if ( fabs( z ) > norma ) norma = fabs( z );
	
	return norma;

}

double wyznacznik( double (&J)[3][3] )
{
	return J[0][0] * J[1][1] * J[2][2] + J[0][1] * J[1][2] * J[2][0] + J[0][2] * J[1][0] * J[2][1]
		- ( J[0][2] * J[1][1] * J[2][0] + J[0][0] * J[1][2] * J[2][1] + J[0][1] * J[1][0] * J[2][2] );
}

void obliczanie_pierwiastków()
{
	double x, y, z, delta_x = 0, delta_y = 0, delta_z = 0, det_x, det_y, det_z, f1, f2, f3;
	double J[3][3], J_pom[3][3], det_J;
	double EST, RESIDIUM;

	x = y = z = 1;

	cout << "  n |             x |             y |             z |            EST |     RESIDIUM |" << endl;
	cout << "----------------------------------------------------------------------------------" << endl;

	for ( int i = 0; i < il_petli; ++i )
	{
		f1 = funkcja_1( x, y, z );
		f2 = funkcja_2( x, y, z );
		f3 = funkcja_3( x, y, z );

		EST = norma( delta_x, delta_y, delta_z );
		RESIDIUM = norma( f1, f2, f3 );

		cout.width( 4 );
		cout << i << "|";
		cout.width( 15 );
		cout << x << "|";
		cout.width( 15 );
		cout << y << "|";
		cout.width( 15 );
		cout << z << "|";
		cout.width( 15 );
		cout << EST << "|";
		cout.width( 15 );
		cout << RESIDIUM << "|" << endl;

		J[0][0] = 2 * x;	J[0][1] = 2 * y;	J[0][2] = 0;
		J[1][0] = 4 * x;	J[1][1] = -1;		J[1][2] = 0;
		J[2][0] = 0.25*(1 + tan( x / 4 )*tan( x / 4 ));
		J[2][1] = 3 * y * y;
		J[2][2] = -5 * z * z * z * z;

		det_J = wyznacznik( J );
	//	cout << "glowny	" << det_J << endl;
		//wyznacnzik x
		J_pom[0][0] = f1;
		J_pom[1][0] = f2;
		J_pom[2][0] = f3;

		for ( int j = 0; j < 3; ++j )
			for ( int k = 1; k < 3; k++ )
				J_pom[j][k] = J[j][k];
			
		det_x = wyznacznik( J_pom );
	//	cout << "x	" << det_x << endl;
		//wyznacnzik y
		J_pom[0][1] = f1;
		J_pom[1][1] = f2;
		J_pom[2][1] = f3;

		for ( int j = 0; j < 3; ++j )
			for ( int k = 0; k < 3; k++ )
				if ( k != 1 )
					J_pom[j][k] = J[j][k];
	
		det_y = wyznacznik( J_pom );
		//cout << "y	" << det_y << endl;
		//wyznacnzik z
		J_pom[0][2] = f1;
		J_pom[1][2] = f2;
		J_pom[2][2] = f3;

		for ( int j = 0; j < 3; ++j )
			for ( int k = 0; k < 2; k++ )
				J_pom[j][k] = J[j][k];

		det_z = wyznacznik( J_pom );
	//	cout << "x	" << det_z << endl;
		delta_x = det_x / det_J;
		delta_y = det_y / det_J;
		delta_z = det_z / det_J;

		x -= delta_x;
		y -= delta_y;
		z -= delta_z;

		
	
		if ( fabs( f1 ) < eps && fabs( f2 ) < eps && fabs( f3 ) < eps ) break;
		
		if ( fabs( x ) < eps && fabs( y ) < eps && fabs( z ) < eps ) break;
		
	}

	
}
int main()
{

	obliczanie_pierwiastków();

	system ("Pause");
	return 0;
}

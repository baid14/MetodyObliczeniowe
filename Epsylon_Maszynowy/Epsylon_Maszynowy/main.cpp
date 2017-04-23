#include <iostream>
#include <iomanip> 
#include <Windows.h>

using namespace std;

void double_epsylon( double liczba )
{
	int t = 0.0;
	double epsylon = liczba / liczba;

	//cout << "Epsylon \t 1+ epsylon \t mantysa." << endl;

	do
	{
		//	cout << epsylon << "\t" << setprecision( 20 ) << 1.0 + epsylon << "\t" << t << endl;
		epsylon *= 0.5;
		t++;

	} while ( ( 1.0 + epsylon/2 ) > 1.0 );	//sprawdzam czy po next przebigu bd mozna dodac do mantysy cos jeszcze

	cout << "Epsylon maszynowy (double):";
	cout.width( 25 );
	cout << epsylon << endl; 
	cout << "Bity przeznaczone na mantyse (double):  " << t << endl << endl;
}

void float_epsylon( float liczba )
{
	int t = 0.0;
	float epsylon = liczba / liczba;
		
	do
	{
		epsylon *= 0.5f;
		t++;

	} while ( ( 1.0f + epsylon / 2 ) > 1.0 );

	cout << "Epsylon maszynowy (float): ";
	cout.width( 25 );
	cout << epsylon << endl;
	cout << "Bity przeznaczone na mantyse (float) :  " << t << endl << endl;
}

int main()
{
	double_epsylon( 20 );
	float_epsylon( 30 );
	
	system( "pause" );

	return 0;
}
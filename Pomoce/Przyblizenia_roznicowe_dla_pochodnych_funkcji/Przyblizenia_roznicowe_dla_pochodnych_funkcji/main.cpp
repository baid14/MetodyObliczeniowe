#define _USE_MATH_DEFINES

#include <iostream>
#include <Windows.h>
#include <cmath>
#include <fstream>
#include <string>
#include <typeinfo>

using namespace std;

const int il_krokow = 20;
const int il_metod = 7;

template <typename T> T roznica_progresywna( T x, T h )
{
	return ( sin( x + h ) - sin( x ) ) / h;
}

template <typename T> T roznica_wsteczna( T x, T h )
{
	return ( sin( x ) - sin( x - h ) ) / h;
}

template <typename T> T roznica_centralna( T x, T h )
{
	return( sin( x + h ) - sin( x - h ) ) / ( 2.0 * h );
}

template <typename T> T roznica_progresywna_trzypunktowa( T x, T h )
{
	return ( 3.0 * sin( x ) + 4.0 * sin( x + h ) - sin( x + 2.0 * h ) ) / ( 2.0 * h );
}

template <typename T> T roznica_wsteczna_trzypunktowa( T x, T h )
{
	return ( sin( x - 2.0 * h ) - 4.0 * sin( x - h ) + 3.0 * sin( x ) ) / ( 2.0 * h );
}

template <typename T> T pochodna( T x )
{
	return cos( x );
}

template <typename T> void oblicz_roznice( T *kroki, T **wyniki )
{
	T poczatek = 0.0;
	T srodek = M_PI_4;
	T koniec = M_PI_2;
	T krok = 0.1;
	//T krok = 0.5;

	for ( int i = 0; i < il_krokow; i++ )
	{
		kroki[i] = krok;

		wyniki[0][i] = fabs( pochodna( poczatek ) - roznica_progresywna( poczatek, krok ) );
		wyniki[1][i] = fabs( pochodna( poczatek ) - roznica_progresywna_trzypunktowa( poczatek, krok ) ); 

		wyniki[2][i] = fabs( pochodna( srodek ) - roznica_progresywna( srodek, krok ) ); 
		wyniki[3][i] = fabs( pochodna( srodek ) - roznica_centralna( srodek, krok ) );
		wyniki[4][i] = fabs( pochodna( srodek ) - roznica_wsteczna( srodek, krok ) ); 

		wyniki[5][i] = fabs( pochodna( koniec ) - roznica_wsteczna( koniec, krok ) ); 
		wyniki[6][i] = fabs( pochodna( koniec ) - roznica_wsteczna_trzypunktowa( koniec, krok ) );

		//krok *= 0.5;
		krok *= 0.1 ;

	}
}

template <typename T> void zapis_danych( T *kroki, T **wyniki, string *opisy )
{
	fstream plik;
	string nazwa = "Dane_";

	nazwa += typeid( T ).name();
	//nazwa += ".dat";
	nazwa += ".txt";
	cout << endl << endl;
	cout << "Zapisuje dane do pliku." << endl;

	plik.open( nazwa.c_str(), fstream::out );

	/*for ( int i = 0; i < il_metod + 1; i++ )
	{
		plik << opisy[i];
		//cout << opisy[i] << "|";
	}*/


	cout << "  n | PP - rpd | PP - rwT | PS - rpd | PS - rcd | PS - rwd | PK - rwd | PK - rwt |"  << endl;
	cout << "----------------------------------------------------------------------------------" << endl;

	plik << endl;

	for ( int i = 0; i < il_krokow; i++ )
	{
		cout.width( 5 );
		//cout << kroki[i]  << "|";
		cout << log10(kroki[i]) << "|";

		plik << log10( kroki[i]) << " ";
		//plik << kroki[i]  << " ";

		for ( int j = 0; j < il_metod ; j++ )
		{
			plik << log10( wyniki[j][i]) << " ";
			//plik <<  wyniki[j][i]  << " ";
			cout.width( 15 );
			//cout <<  wyniki[j][i]  << "|";
			cout << log10(wyniki[j][i]) << "|";
		}
			
		plik << endl;
		cout << endl;
	}

	plik.close();
}

int main()
{

	float *krokiFloat, **wynikiFloat;
	double *krokiDouble, **wynikiDouble;
	string *opisy;

	krokiFloat = new float[il_krokow];
	krokiDouble = new double[il_krokow];

	wynikiFloat = new float *[il_metod];
	wynikiDouble = new double *[il_metod];

	for ( int i = 0; i < il_metod; i++ )
	{
		wynikiFloat[i] = new float [il_krokow];
		wynikiDouble[i] = new double[il_krokow];
	}

	opisy = new string[il_metod + 1];

	/*opisy[0] = "# PP - rpd"; //Punkt poczatkowy - ró¿nica progresywna dwupunktowa
	opisy[1] = "# PP - rwt"; //Punkt poczatkowy - ró¿nica progresywna trzypunktowa

	opisy[2] = "# PS - rpd"; //Punkt œrodkowy - ró¿nica progresywna dwupunktowa
	opisy[3] = "# PS - rcd"; //Punkt œrodkowy - ró¿nica centralna dwupunktowa
	opisy[4] = "# PS - rwd"; //Punkt œrodkowy - ró¿nica wsteczna dwupunktowa

	opisy[5] = "# PK - rwd"; //Punkt koñcowy - ró¿nica wsteczna dwupunktowa
	opisy[6] = "# PK - rwt"; //Punkt koñcowy - ró¿nica wsteczna trzypunktowa*/
	opisy[0] = "Krok,";
	opisy[1] = "PP - progresywna dwupunktowa,"; //
	opisy[2] = "PP - progresywna trzypunktowa,"; //

	opisy[3] = "PS - progresywna dwupunktowa,"; //
	opisy[4] = "PS - centralna dwupunktowa,"; //
	opisy[5] = "PS - wsteczna dwupunktowa,"; //

	opisy[6] = "PK - wsteczna dwupunktowa,"; //
	opisy[7] = "PK - wsteczna trzypunktowa,"; //

	oblicz_roznice( krokiFloat, wynikiFloat );
	oblicz_roznice( krokiDouble, wynikiDouble );

	zapis_danych( krokiFloat, wynikiFloat, opisy );
	zapis_danych( krokiDouble, wynikiDouble, opisy );
	
	system( "Pause" );
	return 0;
}
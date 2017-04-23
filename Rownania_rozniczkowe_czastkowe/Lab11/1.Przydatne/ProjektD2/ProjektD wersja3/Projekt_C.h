//#include "Projekt_C.h"
/*#include <iostream>
//#include <iomanip>
#include <math.h>
#include <Windows.h>
#include <fstream>
//#include <sstream>
#include <string>
#include "calerf.h"


using namespace std;
*/
enum metoda_rozwiazania {THOMAS, LU};

class Projekt_C
{
public:
	Projekt_C(void);
	~Projekt_C(void);
	void informacje(void);										//wy�wietlenie tekstu pocz�tkowego
	void ustaw_zmienne(int n);									//metoda ustawiaj�ca podstawowe zmienne np. krok dla x, krok dla t, lambda
	void zresetowanie();										//metoda zeruj�ca zmienne

	void rozw_anal();											//rozwiazanie analityczne zadanego r�wnania
	void rozw_laasonen(metoda_rozwi�zania metoda);										//rozwi�zanie metod� po�redni� Laasonen 

	double getD() { return D; }
	double getH() { return h; }
	double getR() { return r; }
	double getT_MIN() { return T_MIN; }
	double getT_MAX() { return T_MAX; }
	double getA() { return a; }
	int getN() { return N; }
	int getM() { return M; }
	double getDt() { return dt; }
	double *getPoz() { return pozycja; }
	double *getCzas() { return czas; }
	double **getAnal() { return anal; }
	double **getLaas() { return laasonen; }

private:
	/** Parametry pocz�tkowe */
	const double D;											//wsp�czynnik transportu ciep�a
	const double T_MAX;										//maksymalna czas
	const double T_MIN;										//minimalny czas
	/*************************/

	double a;												//g�rna granica dla wsp�rz�dnej przestrzennej "x"
	double r; 												//poczatek przedzialu zmiennej x
	double h;												//krok wsp�rz�dnej przestrzennej "x"
	double dt;												//krok zmiennej czasu
	double lambda;											//d�zymy do osi�gni�cia LAMBDA = 1 dla metod po�rednich

	int M;													//liczba podzialow przedzialu przestrzennego (zmienna x)				
	int N;													//liczba podzialow przedzialu czasowego (zmienna t)
	double **anal;
	double ** laasonen;

	double *pozycja;
	double *czas;

	void zapisz_do_pliku(char* nazwa_pliku, double** macierz);
	void thomas(double **A, double *b, double *wyn);
	void dekompozycja_LU(double **A, double *b, double *wyn);
};

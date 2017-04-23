#pragma once
#include "stdafx.h"
#include <time.h>

enum metodaRozwi�zania {THOMAS, LU};

class Projekt_D{
	public:
		Projekt_D(void);
		~Projekt_D(void);
		void informacje(void);										//wy�wietlenie tekstu pocz�tkowego
		void ustawZmienne(int n);									//metoda ustawiaj�ca podstawowe zmienne np. krok dla x, krok dla t, lambda
		void zresetowanie();										//metoda zeruj�ca zmienne
		void rozwAnalityczne();										//rozwiazanie analityczne zadanego r�wnania
		void rozwLaasonen(metodaRozwi�zania metoda);				//rozwi�zanie metod� po�redni� Laasonen 
		void rozwKMB(int n);										//rozwi�zanie klasyczn� metod� bezpo�redni�

		// Metody pobieraj�ce dane
		double getD() { return D; }
		double getH() { return h; }
		double getT_MIN() { return T_MIN; }
		double getT_MAX() { return T_MAX; }
		double getA() { return a+1; }
		int getN() { return N; }
		int getM() { return M; }
		double getDt() { return dt; }
		double *getPoz() { return pozycja; }
		double *getCzas() { return czas; }
		double **getAnal() { return anal; }
		double **getLaas() { return laasonen; }
		double **getBezpo(){ return bezpo; }


	private:
		// Parametry pocz�tkowe
		double D;													//wsp�czynnik transportu ciep�a
		const double T_MAX;											//maksymalna czas
		const double T_MIN;											//minimalny czas
		double r;
		double a;													//g�rna granica dla wsp�rz�dnej przestrzennej "x"
		double h;													//krok wsp�rz�dnej przestrzennej "x"
		double dt;													//krok zmiennej czasu
		double lambda;												//d��ymy do osi�gni�cia LAMBDA = 1 dla metod po�rednich
		int M;														//liczba podzialow przedzialu przestrzennego (zmienna x)				
		int N;														//liczba podzialow przedzialu czasowego (zmienna t)
		double **anal;
		double ** laasonen;
		double **bezpo;
		double *pozycja;
		double *czas;

		void zapiszDoPliku(char* nazwaPliku, double** macierz);
		void thomas(double **A, double *b, double *wyn);
		void dekompozycjaLU(double **A, double *b, double *wyn);


};
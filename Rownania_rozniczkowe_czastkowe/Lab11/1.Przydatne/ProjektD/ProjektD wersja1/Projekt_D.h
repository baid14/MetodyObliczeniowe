#pragma once
#include "stdafx.h"
#include <time.h>

enum metodaRozwi¹zania {THOMAS, LU};

class Projekt_D{
	public:
		Projekt_D(void);
		~Projekt_D(void);
		void informacje(void);										//wyœwietlenie tekstu pocz¹tkowego
		void ustawZmienne(int n);									//metoda ustawiaj¹ca podstawowe zmienne np. krok dla x, krok dla t, lambda
		void zresetowanie();										//metoda zeruj¹ca zmienne
		void rozwAnalityczne();										//rozwiazanie analityczne zadanego równania
		void rozwLaasonen(metodaRozwi¹zania metoda);				//rozwi¹zanie metod¹ poœredni¹ Laasonen 
		void rozwKMB(int n);										//rozwi¹zanie klasyczn¹ metod¹ bezpoœredni¹

		// Metody pobieraj¹ce dane
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
		// Parametry pocz¹tkowe
		double D;													//wspó³czynnik transportu ciep³a
		const double T_MAX;											//maksymalna czas
		const double T_MIN;											//minimalny czas
		double r;
		double a;													//górna granica dla wspó³rzêdnej przestrzennej "x"
		double h;													//krok wspó³rzêdnej przestrzennej "x"
		double dt;													//krok zmiennej czasu
		double lambda;												//d¹¿ymy do osi¹gniêcia LAMBDA = 1 dla metod poœrednich
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
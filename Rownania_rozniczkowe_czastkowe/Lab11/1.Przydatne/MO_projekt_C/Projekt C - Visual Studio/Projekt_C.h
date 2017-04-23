#pragma once
#include "stdafx.h"

enum metoda_rozwi¹zania {THOMAS, LU};

class Projekt_C
{
public:
	Projekt_C(void);
	~Projekt_C(void);
	void informacje(void);										//wyœwietlenie tekstu pocz¹tkowego
	void ustaw_zmienne(int n);									//metoda ustawiaj¹ca podstawowe zmienne np. krok dla x, krok dla t, lambda
	void zresetowanie();										//metoda zeruj¹ca zmienne

	void rozw_anal();											//rozwiazanie analityczne zadanego równania
	void rozw_laasonen(metoda_rozwi¹zania metoda);										//rozwi¹zanie metod¹ poœredni¹ Laasonen 

	double getD() { return D; }
	double getH() { return h; }
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
	/** Parametry pocz¹tkowe */
	const double D;											//wspó³czynnik transportu ciep³a
	const double T_MAX;										//maksymalna czas
	const double T_MIN;										//minimalny czas
	/*************************/

	double a;												//górna granica dla wspó³rzêdnej przestrzennej "x"
	double h;												//krok wspó³rzêdnej przestrzennej "x"
	double dt;												//krok zmiennej czasu
	double lambda;											//d¹zymy do osi¹gniêcia LAMBDA = 1 dla metod poœrednich

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
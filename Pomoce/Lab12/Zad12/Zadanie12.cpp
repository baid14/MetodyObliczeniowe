#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
using namespace std;

/*

Metody Obliczeniowe
By Wahuu

*/


double const PI = 3.14;
int const N = 100;	// liczba krokow
int const liczba_wezlow = 7; // liczba wezlow
double const XMIN = -1.0; // poczatek przedzialu
double const XMAX = 1.0; // koniec przedzialu
double const skok_liczba_wezlow = (XMAX - XMIN) / (liczba_wezlow - 1);	//skok dla liczby wezlow
double const skok = (XMAX - XMIN) / N;	//skok dla liczby krokow

// Interpolacja Lagrange
double Lagrange(double x, double *xi, double *yi)
{
    double iloczyn[liczba_wezlow];
    double wynik = 0.0;
    
    for(int i = 0 ; i < liczba_wezlow ; i++)
    {
        iloczyn[i] = 1.0;
        for(int j = 0 ; j < liczba_wezlow ; j++)
            if(j != i)
                iloczyn[i] = iloczyn[i] * ((x - xi[j]) / (xi[i] - xi[j]));
    }
    
    for(int i = 0 ; i < liczba_wezlow ; i++)
        wynik = wynik + yi[i] * iloczyn[i];
        
    return wynik;
}

double Funkcja(double x)
{
	return 1.0 / (1.0 + 25.0 * x * x);
}

double Czebyszew(double x)
{
	double xCz[liczba_wezlow];
    double yCz[liczba_wezlow];
    
    for(int i = 0; i < liczba_wezlow; i++)
    {
        xCz[i] = ((XMAX + XMIN) / 2.0) + ((XMAX - XMIN) / 2.0) * cos(((2.0 * i + 1.0) * PI) / (2.0 * liczba_wezlow)); // wyznaczanie kolejnego wezla
        yCz[i] = Funkcja(xCz[i]); // wyliczanie funkckji dla wezla
    }
	 
    return Lagrange(x, xCz, yCz);
								                         
}

double Rownoodlegle(double x)
{
	double xi[liczba_wezlow];
    double yi[liczba_wezlow];
    
	int i = 0;
	
	for(double z = XMIN; z <= XMAX; z += skok_liczba_wezlow)
    {          
        xi[i] = z; // poczatek przedzialu
        yi[i] = Funkcja(z); // wyliczenie funkcji dla wezla 
        i++;
    }
	return Lagrange(x, xi, yi);
}

int main()
{
    ofstream analitycznie("analitycznie.txt");
    ofstream krok("krok.txt");
    ofstream czebyszew("czebyszew.txt");
    ofstream rownoodlegle("rownoodlegle.txt");
    
	int i = 0;
	double x = XMIN;
	double wynik_Funkcja[N + 1]; 
	double wynik_Czebyszew[N + 1];
	double wynik_Rownoodlegle[N + 1];
    
    cout << "Zadanie 12. \n\n";
    cout << "      krok" << "\t" << "   Analitycznie" << "\t    Czebyszew" << "\t" << "   Rownoodlegle" << endl << endl;
	
	for(int i = 0 ; i < N + 1; i++)
	{			
        wynik_Funkcja[i] = Funkcja(x);
        wynik_Czebyszew[i] = Czebyszew(x);
		wynik_Rownoodlegle[i] = Rownoodlegle(x);
		
		printf("%12lf \t %12lf \t %12lf \t %12lf \n",x,wynik_Funkcja[i],wynik_Czebyszew[i],wynik_Rownoodlegle[i]);
		analitycznie << wynik_Funkcja[i] << endl;
		czebyszew << wynik_Czebyszew[i] << endl;
		rownoodlegle << wynik_Rownoodlegle[i] << endl;
		krok << x << endl;

		x += skok;
	}  
	system("pause");
}
                 

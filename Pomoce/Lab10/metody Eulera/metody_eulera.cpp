#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
using namespace std;

class Rownanie{
public:
	double bezposrednia_metoda_Eulera(double h, double t);
	double posrednia_metoda_Eulera(double h, double t);
	double metoda_trapezow(double h, double t);

	double wzor_analityczny(double t);
};

/**
*	Wyliczenie warto�ci w punkcie ze wzoru analitycznego
*	@return Warto�� dok�adna w punkcie "t" wyliczona ze wzoru analitycznego
*/
double Rownanie::wzor_analityczny(double t)
{
	return 1 - exp( -10.0 * ( t + atan( t ) ));//1+t)*exp(-t);
}

/**
*	Metoda obliczaj�ca przybli�enie r�wnania r�niczkowego pierwszego stopnia w punkcie z danym krokiem bezspo�redni� metod� Eulera
*	@param h Krok do zastosowania w obliczeniach. Im mniejszy tym bardziej dok�adne obliczenia
*	@param t Punkt, dla kt�rego wyliczamy przybli�enie
*	@return Przybli�ona warto�� r�wnania r�zniczkowego w punkcie i z krokiem przekazanym do metody
*/
double Rownanie::bezposrednia_metoda_Eulera(double h, double t)
{
	double blad = 0.0;
	double y = 0.0;			//warunek pocz�tkowy; y(0) = 1
	double wart_dokl = 0.0;

	for(double i=0; i<5.0; i+=h)
	{
		wart_dokl = wzor_analityczny(i);
		//y = y + h * (-y + exp(-i));	//obliczanie kolejnych przybli�e� ze wzoru: x(h+t) = x(t) + h * x'(t)
        y = y + h *( - ((10.0 *i * i + 20.0) /( i * i + 1) )*(y - 1) );
		wart_dokl = fabs(wart_dokl - y);
		if(wart_dokl > blad)
         blad = wart_dokl;
	 }
	return blad;
}

/**
*	Metoda obliczaj�ca przybli�enie r�wnania r�niczkowego pierwszego stopnia w punkcie z danym krokiem po�rednia metod� Eulera
*	@param h Krok do zastosowania w obliczeniach. Im mniejszy tym bardziej dok�adne obliczenia
*	@param t Punkt, dla kt�rego wyliczamy przybli�enie
*	@return Przybli�ona warto�� r�wnania r�zniczkowego w punkcie i z krokiem przekazanym do metody
*/
double Rownanie::posrednia_metoda_Eulera(double h, double t)
{
	double blad = 0.0;
	double y = 1.0;			//warunek pocz�tkowy; y(0) = 1
	double wart_dokl = 0.0;

	for(double i=0; i<t; i+=h)
	{
		wart_dokl = wzor_analityczny(i);
		y = (y + exp(-i) * h) / (1 + h);		//obliczanie kolejnych przybli�e� ze wzoru: x(h+t) = x(t) + h * f(t+1/2h,y(t)+1/2h*y(t, y(t))

		wart_dokl = fabs(wart_dokl - y);
		if(wart_dokl > blad)
         blad = wart_dokl;
	}

	return blad;
}

/**
*	Metoda obliczaj�ca przybli�enie r�wnania r�niczkowego pierwszego stopnia w punkcie z danym krokiem metod� Trapez�w
*	@param h Krok do zastosowania w obliczeniach. Im mniejszy tym bardziej dok�adne obliczenia
*	@param t Punkt, dla kt�rego wyliczamy przybli�enie
*	@return Przybli�ona warto�� r�wnania r�zniczkowego w punkcie i z krokiem przekazanym do metody
*/
double Rownanie::metoda_trapezow(double h, double t)
{
	double blad = 0.0;
	double y = 1.0;			//warunek pocz�tkowy; y(0) = 1
	double wart_dokl = 0.0;

	for(double i=0; i<t; i+=h)
	{
		wart_dokl = wzor_analityczny(i);
		y = (2.0 * y -((y - exp(-i) - exp(-i + h)) * h)) / (2.0 + h);				//obliczanie kolejnych przybli�e� ze wzoru: x(h+t) = x(t) + 1/2 * h * (f(t,y(t))+f(t+h, y(t)+h*f(t,y(t)))

		wart_dokl = fabs(wart_dokl - y);
		if(wart_dokl > blad)
         blad = wart_dokl;
	}

	return blad;
}

int main()
{
	Rownanie rown;
	double krok = 0.1;	//krok dla metod
	const double tmax = 20;		//maks. warto�� punktu
	fstream plik;				//plik do zapisania wynik�w

	plik.open("wyniki.txt", ios::out);
	if( plik.good() != true )
		cout << "Nie uzyskano dostepu do pliku!" << endl;

	plik << "Rozwiazanie dla punktu: t = ;" << tmax << endl << "Wartosc dokladna = ;" << rown.wzor_analityczny(tmax) << endl << "B��dy" << endl;
	plik << "krok;bez. Eulera;pos. Eulera;M. trapezow" << endl;

	for(int i=0; i<20; i++)
	{
		cout <<  krok << ";" << rown.bezposrednia_metoda_Eulera(krok, tmax) <<endl;
			//<<  ";" << rown.posrednia_metoda_Eulera(krok, tmax) << ";" << rown.metoda_trapezow(krok, tmax) << endl;
		krok/=2;
	}

	plik.close();
	return 0;
}

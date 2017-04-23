#include <iostream>
#include <cmath>
using namespace std;
                        //DANIEL KURA GR 23A 98022
class rownanie_kwadratowe{
	double a,b,c;      // pola skladowe klasy odpowiadajace za wartosc wspolczynnikow
	double x0,x1;      //znane miejsca zerowe
	double x0_wyl, x1_wyl, x1_popraw, bw_x0, bw_x1, bw_x1_popraw, delta; // pola sluzace do przechowywania wynikow

public: // gettery poszczegolnych zmiennych
	double getA() { return a; }
	double getB() { return b; }
	double getC() { return c; }
	double getX1() { return x1; }
	double getX0() { return x0; }
	double getX0_wyl() { return x0_wyl; }
	double getX1_wyl() { return x1_wyl; }
	double getX1_popraw() { return x1_popraw; }
	double getBw_x0() { return bw_x0; }
	double getBw_x1() { return bw_x1; }
	double getBw_x1_popraw() { return bw_x1_popraw; }

	rownanie_kwadratowe() {a=1; x1= -0.0295873; x0=2.02959; b=-a*(x0+x1); c=1*x0*x1; }

	void miejsca_zerowe_szkola() // metoda obliczajaca miejsca zerowe metoda "szkolna"
	{
		delta = (b*b) - (4*a*c);
		x0_wyl = (-b + sqrt(delta))/(2*a);
		x1_wyl = (-b - sqrt(delta))/(2*a);
	}

	void miejsce_zerowe_Vieta() // miejsce zerowe obliczanie z wzorow wieta
	{
		x1_popraw = c/(a*x0); // x1*x2 = c/a => x2 = c/(a*x1)
	}

	void bledy_wzgledne() // obliczanie bledow wzglednych
	{
		bw_x0 = (x0 - x0_wyl)/x0;
		bw_x1 = (x1 - x1_wyl)/x1;
		bw_x1_popraw = (x1 - x1_popraw)/x1;
	}
};

int main()
{
	rownanie_kwadratowe przyklad; // stworzenie obiektu

	przyklad.miejsca_zerowe_szkola(); // wykonanie metody oblicznia miejsc zerowerogo "szkola"
	przyklad.miejsce_zerowe_Vieta(); // wykonanie metody oblicznia miejsca zerowerogo "Vieta"
	przyklad.bledy_wzgledne();      // wykonanie metody oblicznia bledow wzglednych

	cout << "Rownanie kwadratowe: " << przyklad.getA() << "x^2 + " << przyklad.getB() << "x + " << przyklad.getC() << endl;
	cout << "Pierwiastki rownania:\nx1: " << przyklad.getX0() <<"\nx2: " << przyklad.getX1() << endl;

	cout << endl << "Pierwiastki wyliczone metoda szkolna:\nx1: " << przyklad.getX0_wyl() << "\nx2: " << przyklad.getX1_wyl() << endl;
	cout << endl << "Bledy wzgledne dla pierwiastkow wyliczonych metoda szkolna:\nx1: " << przyklad.getBw_x0() << "\nx2: " << przyklad.getBw_x1()<< endl;

	cout << endl << "Pierwiastek wyliczony ze wzoru Viete'a (unikciecie redukcji cyfr znaczacych przy odejmowaniu):\nx2: " << przyklad.getX1_popraw() << endl;
	cout << endl << "Blad wzgledny dla pierwiastka wyliczonego ze wzoru Vieta:\nx2: " << przyklad.getBw_x1_popraw() << endl;

	return 0;
}
/* %%%%%%%%%%%%%%%%%% WNIOSKI %%%%%%%%%%%%%%%%%%
Metoda obliczania pierwiastkow rownania tzw. "szkolna"(w ktorej jest odejmowanie)
jest obarczona bledem redukcji, co mozemy zaobserwowac w momencie obliczenia bledu
wzglednego tych pierwiastkow. Obliczony pierwiastek metoda Vieta, w ktorym
nie wystepuje odejmowanie jest zgodny z znanymi wartosciami pierwiastkow, na co
wskazuje blad wzgledny obliczonego pierwiastka.
*/

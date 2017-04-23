#include <iostream>
#include <fstream>
using namespace std;

class Dyskretyzacja{
	int N;						//liczba podprzedzia³ów
	double x1;					//punkty brzegowe
	double x2;
    double h,					//dlugosc przedzia³ow
			*l,					//dolna przek¹tna macierzy trójdiagonalnej
			*d,					//œrodkowa przek¹tna macierzy trójdiagonalnej
			*u,					//górna przek¹tna macierzy trójdiagonalnej	
			*b,					//wektor wyrazow wolnych
			*blad,				//wektor b³edów obliczeñ 
			*x;					//wektor rozwiazañ macierzy

public:
	/* KONSTRUKTORY */
	Dyskretyzacja() { N = 20; x1 =0; x2 = 2; alokacja_pamieci();}
	Dyskretyzacja(int n) { N = n; x1 =0; x2 = 2; alokacja_pamieci();}
	/* DESTRUKTOR */
	~Dyskretyzacja() { dealokacja_pamieci();}

	/* METODY */
	double U(double x);								//wzór analityczny
	double dyskretyzacja_konwencjonalna();	
	double dyskretyzacja_Numerowa();
	void algorytm_Thomasa();
	double najw_blad();
	void ustawN(int wartosc) { N = wartosc; }
	int dajN() { return N; }
	double dajH() { return h; }

	/** Metoda alokuj¹ca pamiêæ pod wektory	*/
	void alokacja_pamieci()
	{
		h = fabs(x2 - x1) / (N - 1);
		l = new double[N]; d = new double[N];
		u = new double[N]; b = new double[N];
		blad = new double[N]; x = new double [N];
	}

	/** Metoad zwalniaj¹ca pamiêc zarezerwowan¹ pod wektory	*/
	void dealokacja_pamieci()
	{
		delete []l; delete []d; delete []u; delete []b; delete []blad; delete []x;
	}

};

/**
*	Funkcja zwracaj¹ca dok³adny wynik
*	@param x Argument funkcji
*	@return Wartoœæ funkcji w punkcie x
*/
double Dyskretyzacja::U(double x)
{
	return (exp(x / 2.0) - exp(2.0 - x / 2.0)) / (1 - exp(2.0));
}

/**
*	Metoda wyszukuj¹ca najwiêkszy b³ad w obliczeniach
*/
double Dyskretyzacja::najw_blad()
{
 double max = 0;
 for(int i = 0; i < N; i++)
	if(fabs(blad[i]) > max)
		max = fabs(blad[i]);
 return max;
}

double Dyskretyzacja::dyskretyzacja_konwencjonalna()
{ 
	 double a = 0.0;

	 dealokacja_pamieci();
	 alokacja_pamieci();

	 l[0] = 0.0;
	 d[0] = 1.0;
	 u[0] = 0.0;
	 b[0] = 1.0;
 
	 for(int i = 1; i < N-1; i++)
	 {
		 l[i] = 1.0 / (h * h);
		 d[i] = -(8.0 + h * h)/(4.0 * h * h);
		 u[i] = 1.0 / (h * h);
		 b[i] = 0.0;
	 }

	 l[N-1] = 0.0;
	 d[N-1] = 1.0;
	 u[N-1] = 0.0;
	 b[N-1] = 0.0;

	 algorytm_Thomasa();
 
	 for(int i=0; i<N; i++)
	 {
		 blad[i] = fabs(x[i] - U(a));
		 a += h;	
	 }
 
	 return najw_blad();
}

double Dyskretyzacja::dyskretyzacja_Numerowa()
{
	 double a = 0.0;

	 dealokacja_pamieci();
	 alokacja_pamieci();

	 l[0] = 0.0;
	 d[0] = 1.0;
	 u[0] = 0.0;
	 b[0] = 1.0;

	 for(int i = 1; i < N - 1; i++)
	 {
		 l[i] = (48.0 - h * h) / (48.0 * h * h);
		 d[i] = (-48.0 - 5.0 * h * h) / (24.0 * h * h);
		 u[i] = (48.0 - h * h) / (48.0 * h * h);
		 b[i] = 0;
	 }

	 l[N-1] = 0.0;
	 u[N-1] = 0.0;
	 d[N-1] = 1.0;
	 b[N-1] = 0.0;

	 algorytm_Thomasa();

	 for(int i = 0; i < N; i++)
	 {
		 blad[i] = fabs(x[i] - U(a));
		 a += h;	
	 }
 
	 return najw_blad();
}

void Dyskretyzacja::algorytm_Thomasa()
{
 double *r = new double[N];       
 
 for(int i = 1; i < N; i++)
     d[i] = d[i] - (l[i] * u[i - 1]) / d[i - 1];
 r[0] = b[0];
 for(int i = 1; i < N; i++)
     r[i] = b[i] - (l[i] * r[i - 1]) / d[i - 1];
 x[N - 1] = r[N - 1] / d[N - 1];
 for(int i = N - 2; i >= 0; i--)
     x[i] = (r[i] - u[i] * x[i + 1]) / d[i];
 
 delete[] r;	
}


int main()
{
	Dyskretyzacja dys;
	fstream plik;

	plik.open("wyniki.xls", ios::out);
	if( plik.good() != true )
		cout << "Nie uzyskano dostepu do pliku!" << endl;

	plik << "Krok;Metoda konwencjonalna;Metoda Numerowa" << endl;

	while(dys.dajN() < 3000)
	{
		plik << dys.dajH() << ";" << dys.dyskretyzacja_konwencjonalna() << ";" << dys.dyskretyzacja_Numerowa() << endl;
		dys.ustawN(dys.dajN()+10);
	}

	plik.close();
	return 0;
}
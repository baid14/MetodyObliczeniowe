#include <iostream>
#include <iomanip>
#include <math.h>
#include <Windows.h>
#include <fstream>
#include <sstream>
#include <string>
#include "calerf.h"


using namespace std;

enum metoda_rozwiazania {THOMAS, LU};

class Projekt_C
{
public:
	Projekt_C(void);
	~Projekt_C(void);
	void informacje(void);										//wyœwietlenie tekstu pocz¹tkowego
	void ustaw_zmienne(int n);									//metoda ustawiaj¹ca podstawowe zmienne np. krok dla x, krok dla t, lambda
	void zresetowanie();										//metoda zeruj¹ca zmienne

	void rozw_anal();											//rozwiazanie analityczne zadanego równania
	void rozw_laasonen(metoda_rozwiazania metoda);				//rozwi¹zanie metod¹ poœredni¹ Laasonen 

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
	/** Parametry pocz¹tkowe */
	const double D;											//wspó³czynnik transportu ciep³a
	const double T_MAX;										//maksymalna czas
	const double T_MIN;										//minimalny czas
	/*************************/

	double a;												//górna granica dla wspó³rzêdnej przestrzennej "x"
	double r; 												//poczatek przedzialu zmiennej x
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


/**
* Konstruktor domyœlny. Ustawienie poczatkowych parametrów i wyzerowanie pozosta³ych sk³adowych.
*/

Projekt_C::Projekt_C() : D(1), T_MAX(2), T_MIN(0)
{
	a = 10;
	h = 0;
	lambda = 0;
	M = 0;
	N = 0;
	r = 1;
}

/** Destruktor */
Projekt_C::~Projekt_C(void)
{
}

/** Metoda uruchamiana na pocz¹tku programu. Wyœwietlenie podstawowych informacji na temat projektu */
void Projekt_C::informacje()
{
	cout << "Projekt C - Metody obliczeniowe" << endl;
	cout << "Dyskretyzacja:\n- Klasyczna metoda bezposrednich \n- Metoda posrednia Lasonnen" << endl;
	cout << "Rozwiazanie algebraicznych ukladow rownan liniowych:\n- Dekompozycja LU macierzy pelnej\n- Algorytm Thomasa" << endl;
	cout << "\nParametry:\n- D = " << getD() << "\n- T_MAX = " << getT_MAX() << endl;
	cout << "\n Przedzial zmiennej przestrzennej x: ["<< getR() <<";" << getA() << "]" << endl;
}

/** Metoda ustawiaj¹ca pozosta³e zmienne potrzebne do rozwi¹zania równania (wyliczenie przedzia³ów i kroków) i alokuje potrzebn¹ pamiêæ
*	@param n Liczba podprzedzia³ów przedzia³u zmiennej przestrzennej "x"
*/
void Projekt_C::ustaw_zmienne(int m)
{
	M = m;
	h = a/(double) (m-1);
	dt = h * h;
	N = (int) ((T_MAX/dt)+1);
	dt = T_MAX / (double)(N-1);

	lambda = D * dt / (h * h);

    anal = new double*[N];
    for (int i=0;i<N;i++)
		anal[i] = new double[M]; 

    laasonen = new double*[N];
    for (int i=0;i<N;i++)
		laasonen[i] = new double[M]; 

	/* Wartoœæ X w danym kroku */
	pozycja = new double[M];
	for(int i = 1; i < M; i++)
		pozycja[i] = i*h;

	/* Wartoœæ czasu T w danym kroku */
	czas = new double[N];
	for(int i = 0; i < N; i++)
		czas[i] = i*dt;

	cout << "Krok h: " << h << "\nKrok t: " << dt << "\nLiczba zbiorow x: " << M << "\nLiczba zbiorow t: " << N << "\n lambda: " << lambda << endl << endl;
}

/** Metoda zeruj¹ce zmienne */
void Projekt_C::zresetowanie()
{
	N = 0;
	M = 0; 
	h = 0;
	dt = 0;
}

/** Metoda rozwi¹zuj¹ca analitycznie równanie rózniczkowe */
void Projekt_C::rozw_anal()
{
     for(int j=0;j<M;j++)
             anal[0][j]=1;							//warunek pocz¹tkowy U(x,0) = 1
             
     for(int i=0;i<N;i++)
             anal[i][M-1]=1;						//warunek brzegowy U(M,t) = 1
             
     for(int i=0;i<N;i++)
             anal[i][0]=0;							//warunek brzegowy U(r,t) = 0

	double t = dt, x = h;
    for (int i=1; i<N; i++)
    {
			for (int j=1; j<M-1; j++)
			{
				
				anal[i][j] = (1-r/x)*erfc((x -r)/ (2*sqrt(D*t))); // rozwiazanie analityczne
				x = x+h;
			}
		x = h;
		t = t+dt;
    }

    zapisz_do_pliku("rozw_analityczne.csv", anal);
}

/** Rozwi¹zanie rownania róczniczkowego poœredni¹ metod¹ Laasonena 
*	@param metoda Metoda do rozwi¹zania uk³adów równañ liniowych. Mo¿e przyjmowaæ "THOMAS" jeœli chcemy wykorzystaæ
*					metode Thomasa albo "LU" jeœli chcemy wykorzystaæ dekompozycje LU macierzy pe³nej.
*/
void Projekt_C::rozw_laasonen(metoda_rozwiazania metoda)
{
	 double x=r,t=0;
     double *b = new double[M];
     double *wyn = new double[M];
	 
	 double **A = new double *[M];
	 for(int i = 0; i < M; i++) A[i] = new double [M];

	 //wyzerowanie macierzy na pocz¹tek
	 for(int i = 0; i < M; i++)
		 for(int j = 0; j < M; j++)
			 A[i][j] = 0;

     for(int j=0;j<M;j++)
          laasonen[0][j]=1;								//warunek pocz¹tkowy U(x,0) = 1

     for( int i = 0; i < N; i++ )
          laasonen[i][0] = 0;							//warunek brzegowy U(M,t) = 1
     for( int i = 0; i < N; i++ ) 
          laasonen[i][M-1] = 1-(r/(r+a))*erfc(a/(2*sqrt(D*t)));							//warunek brzegowy U(0,t) = 0

     /*Wypelnianie macierzy */
     for( int k = 1; k < N; k++ )
     {
		A[0][0] = 1;
		b[0] = 0; //wynika z lewego warunku brzegowego

          for( int i = 1; i < M-1; i++ )
		  {
               A[i][i] = -( 1 + (2*lambda) );
			   A[i][i+1] = lambda;
			   A[i][i-1] = lambda;
			   b[i] = -laasonen[k-1][i];
		  }

		b[M-1] = 1-(r/(r+a))*erfc(a/(2*sqrt(D*t)));	//wynika z prawego warunku brzegowego
		A[M-1][M-1] = 1;

		switch(metoda)	//sprawdzamy jak¹œ metod¹ chcemy wyliczyæ rozwi¹zanie uk³adu równañ
		{
			case THOMAS:
				thomas(A, b, wyn);				//wyliczenie metod¹ Thomasa
				break;
			case LU:
				dekompozycja_LU(A, b, wyn);		//dekompozycja LU
				break; 
		}

        for( int i = 1; i < M-1; i++ )
			laasonen[k][i] = wyn[i];  // zapisanie do kolejnego wiersza wyliczonych wartoœci
     }

	 //zapisanie rozwi¹zañ do pliku
	 switch(metoda)
	{
		case THOMAS:
			zapisz_do_pliku("rozw_laasonen_thomas.csv", laasonen);
			break;
		case LU:
			zapisz_do_pliku("rozw_laasonen_lu.csv", laasonen);
			break;
	}

	//usuwanie zarezerwowanej pamiêci
	for(int i=M-1;i>=0;i--) delete []A[i];
		delete []A;

	delete b, wyn;
}

/** Metoda wykorzstuj¹c¹ metode Thomasa do obliczenia rozwi¹zañ uk³adu równañ liniowych
*	@param A Macierz wspó³czynników
*	@param b Wektor wyrazów wolnych
*	@param wyn Wektor rozwi¹zañ uk³adu równañ */
void Projekt_C::thomas(double **A, double *b, double *wyn)
{
	double *beta, *gamma;
	beta = new double[M]; gamma = new double[M];

	//obliczanie gammy i bety
	beta[0] = - (A[0][1]/A[0][0]);
	gamma[0] = (b[0]/A[0][0]);

	for(int i=1; i<M; i++)
	{
		if (i<=M-2)
			beta[i] = - ((A[i][2+i-1])/((A[i][i-1]*beta[i-1])+A[i][1+i-1]));

		gamma[i] = (b[i] - (A[i][i-1] * gamma[i-1]))/(A[i][i-1] * beta[i-1] + A[i][1+i-1]);
	}
	beta[M-1] = 0;

	//obliczanie rozwi¹zania
	wyn[M-1]=gamma[M-1];
	for(int i = M - 2; i>=0; i--)
		wyn[i] = beta[i]*wyn[i+1]+gamma[i];

	//zwolnienie zarezerwowanej pamiêci
	delete beta, gamma;
}


/** Metoda wykorzystuj¹c¹ dekompozycje LU metod¹ elimiacji Gaussa do wyliczenia rozwi¹zañ uk³adu równañ
*	@param A Macierz do dekompozycji
*	@param b Wektor wyrazów wolnych
*	@param wyn Wektor rozwi¹zañ uk³adu równañ */
void Projekt_C::dekompozycja_LU(double **A, double *b, double *wyn)
{
	//Dekompozycja LU macierzy- metoda eliminacji Gaussa 
	double x;
	for(int k=0; k<M-1; k++)
	{
		for(int i=k+1; i<M; i++)
		{
			x = A[i][k]/A[k][k];
			A[i][k] = x;
			for(int j=k+1; j<M; j++)
			{
				A[i][j] = A[i][j] - (x*A[k][j]);
			}
		}
	}

	//Rozwi¹zywanie uk³adu równañ
	double suma;
	double *z = new double[M];

	//podstawianie w przód
	for(int i=0; i<M; i++)
	{
		suma = 0;
		for(int j=0; j<=i-1; j++)
		{
			suma += A[i][j]*z[j];
		}

		z[i] = b[i]-suma;
	}

	 //podstawianie w ty³
	for(int i=M-1; i>=0; i--)
	{
		suma = 0;
		for(int j=i+1; j<M; j++)
		{
			suma +=A[i][j]*wyn[j];
		}

		wyn[i] = (z[i]-suma)/A[i][i];
	}
}

/**
*	Metoda zapisuj¹ca macierz wynikow¹ do pliku
*	@param nazwa_pliku Nazwa pod jak¹ chcemy zapisaæ
*	@param macierz Macierz wynikow¹, któr¹ chcemy zapisaæ do pliku
*/
void Projekt_C::zapisz_do_pliku(char *nazwa_pliku, double** macierz)
{
	fstream plik;
	plik.open(nazwa_pliku, ios::out);

	if( plik.good() == true )
	{
		plik << ";";

		for(int i = 0; i<M; i++)
			plik << pozycja[i]<<";";

		plik << endl;

		for(int i=0;i<N;i++)
		{
			plik << czas[i]<<";";
			for(int j=0;j<M;j++)
			{
					plik<<macierz[i][j]<<";";
			}
			plik<<"\n";
		}
	} else cout << "Nie uzyskano dostepu do pliku " << nazwa_pliku << endl;
	
	plik.close();
}




//void wykres_porownujacy(Projekt_C proj, double czas, int czestotliwosc, int metoda);
//void wykres_max_blad(Projekt_C proj, int freq);

int main()
{
	Projekt_C proj;
	proj.informacje();

	proj.ustaw_zmienne(60);
	proj.rozw_anal();
	//proj.rozw_laasonen(THOMAS);
	
	proj.rozw_laasonen(LU);
	system("PAUSE");
	return 0;
	//wykres_porownujacy(proj, proj.getDt(), 20, 1);
	//wykres_porownujacy(proj, proj.getDt(), 20, 2);

	//wykres_max_blad(proj, 50);
//}

/*void wykres_porownujacy(Projekt_C proj, double czas, int czestotliwosc, int metoda)
{
	if ( czas >= proj.getT_MIN() && czas <= proj.getT_MAX())
	{
			int j = (int) (czas * proj.getN());

			double dt = proj.getDt();
			//czas = j*dt;
			int M = proj.getM();
			double *pozycja;
			pozycja = proj.getPoz();

			ostringstream czas2;
			//czas2 << czas;
			
			ostringstream czestotliwosc2;
			czestotliwosc2 << czestotliwosc;
			
			string nazwa_pliku;
			double **x;
			if (metoda == 1)
			{
				x = proj.getAnal();
				nazwa_pliku = "wart_anal ";
			}
			else if (metoda == 2)
			{
				x = proj.getLaas();
				nazwa_pliku = "wart_lassonen ";
			}

			//nazwa_pliku += "-" + czas2.str();
			nazwa_pliku += "-" + czestotliwosc2.str() + ".csv";

			
			double *ptrX;
			ptrX = x[j];

			fstream dane_wykres;
			dane_wykres.open(nazwa_pliku.c_str(), ios::out);
			if ( dane_wykres.is_open() )
			{
				for(int i = 0; i < M; i+= czestotliwosc)
					dane_wykres << pozycja[i] << "\t"  << ptrX[i] << endl; 
				dane_wykres.flush();
				dane_wykres.close();
			}
			else
				cout << endl << "Nie udalo sie stworzyc pliku do zapisu: " << nazwa_pliku;
	}
	else
		cout << endl << "Czas musi zawierac sie w przedziale " << proj.getT_MIN() << " ; " << proj.getT_MAX() << endl;
}

void wykres_max_blad(Projekt_C proj, int freq)
{
		double dt = proj.getDt();
		int N = proj.getN();
		int M = proj.getM();
		double *pozycje;
		pozycje = proj.getPoz();
		double *czasy;
		czasy = proj.getCzas();

		fstream wyniki;
		wyniki.open("wykres3.csv", ios::out);

		if(wyniki.is_open())
		{
			double **x;
			x = proj.getLaas();
			double *ptrX;

			double **an;
			an = proj.getAnal();

			double err;
			double max;

			for( int j = 0; j < N; j+= freq)
			{
				max = 0;
				ptrX = x[j];
				for(int i = 0; i < M; i++)
				{
					err = ptrX[i] - an[j][i];
					err = fabs( err );
					if ( err > max )
						max = err;
				} 
				wyniki << czasy[j] << ";" << max << endl; 
			}

			wyniki.close();
		}
*/
}

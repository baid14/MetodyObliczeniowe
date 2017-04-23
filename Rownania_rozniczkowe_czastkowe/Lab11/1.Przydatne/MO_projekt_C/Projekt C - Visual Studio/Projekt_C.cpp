#include "Projekt_C.h"


/**
* Konstruktor domy�lny. Ustawienie poczatkowych parametr�w i wyzerowanie pozosta�ych sk�adowych.
*/
Projekt_C::Projekt_C() : D(1), T_MAX(2), T_MIN(0)
{
	a = 6 * sqrt(D*T_MAX);
	h = 0;
	lambda = 0;
	M = 0;
	N = 0;
}

/** Destruktor */
Projekt_C::~Projekt_C(void)
{
}

/** Metoda uruchamiana na pocz�tku programu. Wy�wietlenie podstawowych informacji na temat projektu */
void Projekt_C::informacje()
{
	cout << "Projekt C - Metody obliczeniowe" << endl;
	cout << "Dyskretyzacja:\n- Metoda posrednia Lasonnen" << endl;
	cout << "Rozwiazanie algebraicznych ukladow rownan liniowych:\n- Dekompozycja LU macierzy pelnej\n- Algorytm Thomasa" << endl;
	cout << "\nParametry:\n- D = " << getD() << "\n- T_MAX = " << getT_MAX() << endl;
	cout << "\n Przedzial zmiennej przestrzennej x: [0;" << getA() << "]" << endl;
}

/** Metoda ustawiaj�ca pozosta�e zmienne potrzebne do rozwi�zania r�wnania (wyliczenie przedzia��w i krok�w) i alokuje potrzebn� pami��
*	@param n Liczba podprzedzia��w przedzia�u zmiennej przestrzennej "x"
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

	/* Warto�� X w danym kroku */
	pozycja = new double[M];
	for(int i = 0; i < M; i++)
		pozycja[i] = i*h;

	/* Warto�� czasu T w danym kroku */
	czas = new double[N];
	for(int i = 0; i < N; i++)
		czas[i] = i*dt;

	cout << "Krok h: " << h << "\nKrok t: " << dt << "\nLiczba zbiorow x: " << M << "\nLiczba zbiorow t: " << N << "\n lampda: " << lambda << endl << endl;
}

/** Metoda zeruj�ce zmienne */
void Projekt_C::zresetowanie()
{
	N = 0;
	M = 0; 
	h = 0;
	dt = 0;
}

/** Metoda rozwi�zuj�ca analitycznie r�wnanie r�zniczkowe */
void Projekt_C::rozw_anal()
{
     for(int j=0;j<M;j++)
             anal[0][j]=0;							//warunek pocz�tkowy U(x,0) = 0
             
     for(int i=0;i<N;i++)
             anal[i][M-1]=0;						//warunek brzegowy U(M,t) = 0
             
     for(int i=0;i<N;i++)
             anal[i][0]=1;							//warunek brzegowy U(0,t) = 1

	double t = dt, x = h;
    for (int i=1; i<N; i++)
    {
			for (int j=1; j<M-1; j++)
			{
				anal[i][j] = erfc(x / (2*sqrt(D*t))); // rozwiazanie analityczne
				x = x+h;
			}
		x = h;
		t = t+dt;
    }

    zapisz_do_pliku("rozw_analityczne.csv", anal);
}

/** Rozwi�zanie rownania r�czniczkowego po�redni� metod� Laasonena 
*	@param metoda Metoda do rozwi�zania uk�ad�w r�wna� liniowych. Mo�e przyjmowa� "THOMAS" je�li chcemy wykorzysta�
*					metode Thomasa albo "LU" je�li chcemy wykorzysta� dekompozycje LU macierzy pe�nej.
*/
void Projekt_C::rozw_laasonen(metoda_rozwi�zania metoda)
{
	 double x=0,t=0;
     double *b = new double[M];
     double *wyn = new double[M];
	 
	 double **A = new double *[M];
	 for(int i = 0; i < M; i++) A[i] = new double [M];

	 //wyzerowanie macierzy na pocz�tek
	 for(int i = 0; i < M; i++)
		 for(int j = 0; j < M; j++)
			 A[i][j] = 0;

     for(int j=0;j<M;j++)
          laasonen[0][j]=0;								//warunek pocz�tkowy U(x,0) = 0

     for( int i = 0; i < N; i++ )
          laasonen[i][0] = 1;							//warunek brzegowy U(M,t) = 0
     for( int i = 0; i < N; i++ ) 
          laasonen[i][M-1] = 0;							//warunek brzegowy U(0,t) = 1

     /*Wypelnianie macierzy */
     for( int k = 1; k < N; k++ )
     {
		A[0][0] = 1;
		b[0] = 1; //wynika z lewego warunku brzegowego

          for( int i = 1; i < M-1; i++ )
		  {
               A[i][i] = -( 1 + (2*lambda) );
			   A[i][i+1] = lambda;
			   A[i][i-1] = lambda;
			   b[i] = -laasonen[k-1][i];
		  }

		b[M-1] = 0;	//wynika z prawego warunku brzegowego
		A[M-1][M-1] = 1;

		switch(metoda)	//sprawdzamy jak�� metod� chcemy wyliczy� rozwi�zanie uk�adu r�wna�
		{
			case THOMAS:
				thomas(A, b, wyn);				//wyliczenie metod� Thomasa
				break;
			case LU:
				dekompozycja_LU(A, b, wyn);		//dekompozycja LU
				break; 
		}

        for( int i = 1; i < M-1; i++ )
			laasonen[k][i] = wyn[i];  // zapisanie do kolejnego wiersza wyliczonych warto�ci
     }

	 //zapisanie rozwi�za� do pliku
	 switch(metoda)
	{
		case THOMAS:
			zapisz_do_pliku("rozw_laasonen_thomas.csv", laasonen);
			break;
		case LU:
			zapisz_do_pliku("rozw_laasonen_lu.csv", laasonen);
			break;
	}

	//usuwanie zarezerwowanej pami�ci
	for(int i=M-1;i>=0;i--) delete []A[i];
		delete []A;

	delete b, wyn;
}

/** Metoda wykorzstuj�c� metode Thomasa do obliczenia rozwi�za� uk�adu r�wna� liniowych
*	@param A Macierz wsp�czynnik�w
*	@param b Wektor wyraz�w wolnych
*	@param wyn Wektor rozwi�za� uk�adu r�wna� */
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

	//obliczanie rozwi�zania
	wyn[M-1]=gamma[M-1];
	for(int i = M - 2; i>=0; i--)
		wyn[i] = beta[i]*wyn[i+1]+gamma[i];

	//zwolnienie zarezerwowanej pami�ci
	delete beta, gamma;
}


/** Metoda wykorzystuj�c� dekompozycje LU metod� elimiacji Gaussa do wyliczenia rozwi�za� uk�adu r�wna�
*	@param A Macierz do dekompozycji
*	@param b Wektor wyraz�w wolnych
*	@param wyn Wektor rozwi�za� uk�adu r�wna� */
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

	//Rozwi�zywanie uk�adu r�wna�
	double suma;
	double *z = new double[M];

	//podstawianie w prz�d
	for(int i=0; i<M; i++)
	{
		suma = 0;
		for(int j=0; j<=i-1; j++)
		{
			suma += A[i][j]*z[j];
		}

		z[i] = b[i]-suma;
	}

	 //podstawianie w ty�
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
*	Metoda zapisuj�ca macierz wynikow� do pliku
*	@param nazwa_pliku Nazwa pod jak� chcemy zapisa�
*	@param macierz Macierz wynikow�, kt�r� chcemy zapisa� do pliku
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
#pragma once
#include "stdafx.h"
#include "Projekt_D.h"



// Konstruktor domyœlny. Ustawienie poczatkowych parametrów i wyzerowanie pozosta³ych sk³adowych.
Projekt_D::Projekt_D() : T_MAX(2), T_MIN(0){
	a = 10;
	h = 0;															// DO SPRAWDZENIA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	lambda = 0;
	M = 0;
	N = 0;
	r=1;
	D=1;
}

// Destruktor
Projekt_D::~Projekt_D(void){}

// Funkcja wyœwietlaj¹ca podstawowe informacje na temat projektu
void Projekt_D::informacje(){
	cout << "Projekt D - Metody obliczeniowe" << endl;
	cout << "Dyskretyzacja:\n- Metoda posrednia Lasonnen" << endl;
	cout << "Rozwiazanie algebraicznych ukladow rownan liniowych:\n- Dekompozycja LU macierzy pelnej\n- Algorytm Thomasa" << endl;
	cout << "\nParametry:\n- D = " << getD() << "\n- T_MAX = " << getT_MAX() << endl;
	cout << "\n Przedzial zmiennej przestrzennej x: [1;" << getA() << "]" << endl;
}

// Funkcja ustawiaj¹ca pozosta³e zmienne potrzebne do rozwi¹zania równania
void Projekt_D::ustawZmienne(int m){
	M = m;
	h = a/(double) (m-1);
	dt = h * h;															//DO SPRAWDZENIA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	N = (int) ((T_MAX/dt)+1);
	dt = T_MAX / (double)(N-1);
	lambda = D * dt / (h * h);

    anal = new double*[N];												// deklaracja macierzy z rozwi¹zaniami analitycznymi i alokacja pamiêci
    for (int i=0;i<N;i++)
		anal[i] = new double[M]; 

    laasonen = new double*[N];											// deklaracja macierzy z rozwi¹zaniami z metody Laasonen i alokacja pamiêci	
    for (int i=0;i<N;i++)
		laasonen[i] = new double[M];

	for(int i=0;i<N;i++){												// pocz¹tkowe wype³nienie obu macierzy zerami
		for(int j=0;j<M;j++){
			anal[i][j]=0;
			laasonen[i][j]=0;
		}
	}
				
	pozycja = new double[M];											// alokacja pamiêci i wype³nienie tablicy z wartoœciami x w danym kroku
	for(int i = 0; i < M+1; i++)
		pozycja[i] = (i)*h+1;
	
	czas = new double[N];												// alokacja pamiêci i wype³nienie tablicy z wartoœciami t w danym kroku
	for(int i = 0; i < N+1; i++)
		czas[i] = i*dt;

	cout << "Krok h: " << h << "\nKrok t: " << dt << "\nLiczba zbiorow x: " << M << "\nLiczba zbiorow t: " << N << "\n lampda: " << lambda << endl << endl;
}

//Funkcja zeruj¹ca zmienne
void Projekt_D::zresetowanie(){
	N = 0;
	M = 0; 
	h = 0;
	dt = 0;
}

// Funkcja rozwi¹zuj¹ca równanie ró¿niczkowe analitycznie
void Projekt_D::rozwAnalityczne(){
	double t = dt;
	double x = h+1;
     
	for(int j=0;j<M;j++)
		anal[0][j]=1;													// warunek pocz¹tkowy U(x,0) = 1
             
    for(int i=0;i<N;i++)
        anal[i][M-1]=1.0-(r/(r+a))*erfc(a/(2.0*sqrt(D*t)));				// warunek brzegowy U(M,t) = funkcja
             
    for(int i=0;i<N;i++)
        anal[i][0]=0;													// warunek brzegowy U(0,t) = 0

    for (int i=1; i<N; i++){
		for (int j=1; j<M-1; j++){
			anal[i][j] = 1.0-(r/x)*erfc((x-r)/(2.0*sqrt(D*t)));			// rozwiazanie analityczne
			x = x+h;
		}
		x = h+1;
		t = t+dt;
    }

    zapiszDoPliku("rozw_analityczne.csv", anal);
}

// Funkcja rozwi¹zuj¹ca równianie róczniczkowe poœredni¹ metod¹ Laasonena 
// zmiennna 'metoda' - metoda do rozwi¹zania uk³adów równañ liniowych. 
// Mo¿e przyjmowaæ "THOMAS" jeœli chcemy wykorzystaæ metode Thomasa, 
// LU" jeœli chcemy wykorzystaæ dekompozycje LU macierzy pe³nej.
void Projekt_D::rozwLaasonen(metodaRozwi¹zania metoda){
	 double x=1,t=0;
     double *b = new double[M];
     double *wyn = new double[M];

	 double **A = new double *[M];
	 for(int i = 0; i < M; i++) 
		 A[i] = new double [M];

				
	 for(int i = 0; i < M; i++)											// pocz¹tkowe wyzerowanie macierzy
		 for(int j = 0; j < M; j++)
			 A[i][j] = 0;

     for(int j=0;j<M;j++)
          laasonen[0][j]=1;												// warunek pocz¹tkowy U(x,0) = 0

     for( int i = 0; i < N; i++ )
          laasonen[i][0] = 0;											// warunek brzegowy U(M,t) = 0
    
	 for( int i = 0; i < N; i++ ) 
          laasonen[i][M-1] = 1.0-(r/(r+a))*erfc(a/(2.0*sqrt(D*t)));		// warunek brzegowy U(0,t) = 1

     
     for( int k = 1; k < N; k++ )										// wypelnianie macierzy wynikami
     {
		A[0][0] = 1;
		b[0] = 0;														// wynika z lewego warunku brzegowego
        for( int i = 1; i < M-1; i++ ){
			A[i][i] = -( 1 + (2*lambda) );
			A[i][i+1] = lambda;
			A[i][i-1] = lambda;
			b[i] = -laasonen[k-1][i];
		}

		b[M-1] = 1.0-(r/(r+a))*erfc(a/(2.0*sqrt(D*t)));					// wynika z prawego warunku brzegowego
		A[M-1][M-1] = 1;

		switch(metoda){													// wybór w³aœciwej metody do obliczeñ ( Thomas, LU)
			case THOMAS:
				thomas(A, b, wyn);										// obliczenia wykonywane metod¹ Thomasa
				break;
			case LU:
				dekompozycjaLU(A, b, wyn);								// obliczenia wykonywane metod¹ dekompozycji LU
				break; 
		}

        for( int i = 1; i < M-1; i++ )
			laasonen[k][i] = wyn[i];									// zapisanie do kolejnego wiersza wyliczonych wartoœci
     }

	 switch(metoda){													// zapisanie rozwi¹zañ do pliku dla odpowiednich metod (Thomas, LU)
		case THOMAS:
			zapiszDoPliku("rozw_laasonen_thomas.csv", laasonen);
			break;
		case LU:
			zapiszDoPliku("rozw_laasonen_lu.csv", laasonen);
			break;
	}


	for(int i=M-1;i>=0;i--)												// zwalnianie zalokowanej pamiêci
		delete []A[i];
		
	delete []A;
	delete b, wyn;
}

// Funkcja obliczaj¹ca wartoœci metod¹ Thomasa
// A - macierz wspó³czynników
// b - wektor wyrazów wolnych
// wyn - wektor rozwi¹zañ uk³adu równañ
void Projekt_D::thomas(double **A, double *b, double *wyn){
	double *beta;
	double *gamma;
	beta = new double[M]; 
	gamma = new double[M];


	beta[0] = - (A[0][1]/A[0][0]); 										// obliczanie bety
	gamma[0] = (b[0]/A[0][0]);											// obliczanie gammy

	for(int i=1; i<M; i++){
		if (i<=M-2)
			beta[i] = - ((A[i][2+i-1])/((A[i][i-1]*beta[i-1])+A[i][1+i-1]));

		gamma[i] = (b[i] - (A[i][i-1] * gamma[i-1]))/(A[i][i-1] * beta[i-1] + A[i][1+i-1]);
	}

	beta[M-1] = 0;


	wyn[M-1]=gamma[M-1];												// obliczanie rozwi¹zania
	for(int i = M - 2; i>=0; i--)
		wyn[i] = beta[i]*wyn[i+1]+gamma[i];

	
	delete beta, gamma;													// zwalnianie zalokowanej pamiêci
}


// Funkcja obliczaj¹ca wartoœci metod¹ dekompozycji LU
// A - macierz do dekompozycji
// b - wektor wyrazów wolnych
// wyn - wektor rozwi¹zañ uk³adu równañ
void Projekt_D::dekompozycjaLU(double **A, double *b, double *wyn){
	double x;

	for(int k=0; k<M-1; k++){											// metoda eliminacji Gaussa	
		for(int i=k+1; i<M; i++){
			x = A[i][k]/A[k][k];
			A[i][k] = x;
			for(int j=k+1; j<M; j++){
				A[i][j] = A[i][j] - (x*A[k][j]);
			}
		}
	}

	double suma;
	double *z = new double[M];


	for(int i=0; i<M-1; i++){											//podstawianie w przód	
		for(int j=i+1;j<M;j++){
			b[j]=b[j]-A[j][i]*b[i];
		}
	}

	 
	for(int i=M-1; i>=0; i--){											//podstawianie w ty³	
		suma = 0;
		for(int j=i+1; j<M; j++){
			suma +=A[i][j]*wyn[j];
		}
		wyn[i] = (b[i]-suma)/A[i][i];
	}
}


//	Funkcja zapisuj¹ca macierz do pliku
void Projekt_D::zapiszDoPliku(char *nazwaPliku, double** macierz){
	fstream plik;
	plik.open(nazwaPliku, ios::out);

	if( plik.good() == true ){
		plik << ";";

		for(int i = 0; i<M; i++)
			plik << pozycja[i]<<";";

		plik << endl;

		for(int i=0;i<N;i++){
			plik << czas[i]<<";";
			for(int j=0;j<M;j++){
					plik<<macierz[i][j]<<";";
			}
			plik<<"\n";
		}
	} else 
		cout << "Nie uzyskano dostepu do pliku " << nazwaPliku << endl;
	
	plik.close();
}


void Projekt_D::rozwKMB(int n){
	lambda = 0.4;
	double  t = 0,x = 1;
	double modR = r/h;

	bezpo = new double*[N];
	for (int i = 0; i < N; i++)
		bezpo[i] = new double[M];

	for(int j = 0; j<M; j++){											// warunek pocz¹tkowy
		bezpo[0][j]=1;
	}

	for( int i = 0; i < N; i++ ){										// warunek brzegowy dla z pocz¹tkowego
		bezpo[i][0] = 0;
	}

	for( int i = 0; i < N; i++ ){										// warunek brzegowy dla x koñcowego
		bezpo[i][n-1] = 1.0- (r/(r+a) * erfc(a/(2.0*sqrt(D*t))));
}

double *ptrX;

	ptrX = bezpo[0];
	for( int j = 1; j < N; j++ ) {
		for(int i = 1; i < M-1; i++)
			bezpo[j][i] = (lambda-lambda/(i+modR))*ptrX[i-1] + (1.0-2.0*lambda)*ptrX[i] + (lambda+lambda/(i+modR))*ptrX[i+1];
				
	ptrX = bezpo[j];
	}
zapiszDoPliku ("rozw_klasyczna_met_bezposrednia.csv", bezpo);
}
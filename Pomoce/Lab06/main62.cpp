#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#define SIZE 6

using namespace std;

/* Wyswietlanie macierzy */
void pokaz_macierz(double aa[SIZE][SIZE]){
	for(int i=0;i<SIZE;i++){
		for(int j=0;j<SIZE;j++){
			cout<<setw(10)<<aa[i][j];
		}
        cout<<endl;
	}
}
/* Funkcja wykonujaca operacje na macierzy */
void thomas(double A[SIZE][SIZE]){
    for(int i=1;i<SIZE;i++){
        A[i][i]=A[i][i]-A[i][i-1]*(A[i-1][i]/A[i-1][i-1]);
        A[i][i-1]=A[i][i-1]/A[i-1][i-1];
    }
}
/* Funkcja rozwiazujaca uklad */
void rozwiazanie(double A[SIZE][SIZE],double *b,double *x){
	for(int i=1;i<SIZE;i++){
		b[i]=b[i]-A[i][i-1]*b[i-1];
	}

	for(int i=0;i<SIZE;i++){
			 cout<<" "<<b[i];
	}
	x[SIZE-1]=b[SIZE-1]/A[SIZE-1][SIZE-1];
	for(int i=SIZE-2;i>=0;i--){
    x[i]=(b[i]-A[i][i+1]*x[i+1])/A[i][i];
	}
}

int main(){
	//zmienne pomocnicze
	int i, j;
	//macierz wejsciowa
	double A[SIZE][SIZE], b[SIZE];
	A[0][0] = 1.0; A[0][1] =  1.0 / 2.0; A[0][2]=  0.0; A[0][3] =  0.0; A[0][4] =  0.0; A[0][5] =  0.0;
	A[1][0] =  2.0 / 3.0; A[1][1] = 2.0; A[1][2]=  3.0 / 4.0; A[1][3] =  0.0; A[1][4] =  0.0; A[1][5] =  0.0;
	A[2][0] =  0.0; A[2][1] =  4.0 / 5.0; A[2][2]= 3.0; A[2][3] =  5.0 / 6.0; A[2][4] =  0.0;  A[2][5] =  0.0;
	A[3][0] =  0.0; A[3][1] =  0.0; A[3][2]=  6.0 / 7.0; A[3][3] = 4.0; A[3][4] =  7.0 / 8.0; A[3][5] =  0.0;
	A[4][0] =  0.0; A[4][1] =  0.0; A[4][2]=  0.0; A[4][3] =  8.0 / 9.0; A[4][4] = 5.0; A[4][5] =  9.0 / 10.0;
	A[5][0] =  0.0; A[5][1] =  0.0; A[5][2]=  0.0; A[5][3] =  0.0; A[5][4] = 10.0 / 11.0;; A[5][5] =  6.0;


	b[0] = 7.0 / 2.0; b[1] = -3.0; b[2] = 11.0 / 2.0; b[3] = -191.0 / 28.0; b[4] = 193.0 / 30.0; b[5] = -46.0 / 11.0;


	//wektor b

	//wektor rozwiazania
	double x[SIZE];

	cout<<"Macierz A:"<<endl;
	pokaz_macierz(A);
	cout<<endl<<"----------------------------------------------------"<<endl;
	cout<<"Wektor b:"<< endl;
	for(i=0; i<SIZE; i++){
			 cout <<" "<<b[i];
	}
	cout<<endl<<"----------------------------------------------------"<<endl;
	thomas(A);
	cout<<endl<<"Zastosowanie alg. Thomasa:"<<endl;
	pokaz_macierz(A);
	cout<<endl<<"----------------------------------------------------"<<endl;
	rozwiazanie(A, b, x);
	//wypisanie wyniku
	cout<<endl<<"Ostateczny wynik:"<<endl;
	for(i=0;i<SIZE;i++){
			 cout<<" "<<x[i];
	}
	cout<<endl<<"----------------------------------------------------"<<endl;
	cout<<endl;
	system("pause");
	return 0;
}

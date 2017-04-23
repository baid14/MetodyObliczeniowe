#include <iostream>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

/*

Metody Obliczeniowe
By Wahuu

*/


/* WYSWIETLANIE MACIERZY */
void wyswietl(double** mac, int n, int m){
	for(int i=0;i<n;i++){
		printf("\n|");
		for(int j=0;j<m;j++)
			printf(" %.1lf  ",mac[i][j]);
		printf("|");
	}
}

/* ALGORYTM THOMASA */
double *thomas(double **mac, double *b, int n){
	double *beta, *gamma;
	beta = new double[n]; 
	gamma = new double[n];
	double *x = new double[n];
	//obliczanie bety i gammy
	beta[0] = -(mac[0][1] / mac[0][0]);
	gamma[0] = (b[0] / mac[0][0]);

	for(int i=1;i<n;i++){
		if (i<=n-2)
			beta[i] = - ((mac[i][i+1]) / ((mac[i][i-1] * beta[i-1]) + mac[i][i]));
		gamma[i] = (b[i] - (mac[i][i-1] * gamma[i-1]))/(mac[i][i-1] * beta[i-1] + mac[i][i]);
	}
	beta[n-1] = 0;
	//obliczanie x
	x[n-1]=gamma[n-1];
	for(int i=3;i>=0;i--)
		x[i] = beta[i] * x[i+1] + gamma[i];
	return x;
}

/* PROGRAM G£ÓWNY */
int main(){
	double b[5] = {18.0,32.0,44.0,36.0,22.0};
    double A[5][5] = {10.0,4.0,0.0,0.0,0.0,
                      1.0,11.0,3.0,0.0,0.0,
                      0.0,2.0,12.0,2.0,0.0,
                      0.0,0.0,3.0,13.0,1.0,
                      0.0,0.0,0.0,4.0,14.0};

	double** tabA = new double *[5];
	for (int i=0;i<5;i++) {
		tabA[i] = new double[5];
		memcpy(tabA[i],A[i],sizeof(double)*5);
	}
            
	double* tabB = new double [5];           
	memcpy(tabB,b,sizeof(double)*5);

	double* x = new double [5];

	printf("Macierz troj-diagonalna: ");
	wyswietl(tabA,5,5);

	printf("\nWektor b\n");
	for(int i=0;i<5;i++)
		printf("| %.1lf |\n",tabB[i]);

	x = thomas(tabA, tabB, 5);

	printf("\nRozwiazanie: \n");
	printf("Wektor x\n");
	for(int i=0;i<5;i++)
		printf("| %.1lf |\n",x[i]);
		system("pause");
}

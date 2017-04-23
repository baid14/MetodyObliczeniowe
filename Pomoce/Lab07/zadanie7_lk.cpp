#include <iostream>
#include <math.h>
#include <iomanip>
#include <stdlib.h>

using namespace std;

const int dimM = 3;
const int dimN = 3;
const double eps = 0.0000001;
const int loop = 50;

void showMatrix(double a[][dimN]) {
	for(int i = 0; i < dimM; i++) {
		for(int j = 0; j < dimN; j++)
			cout << " " << a[i][j];
             cout << endl;
	}
}

void showVector(double v[dimM]) {
	for(int i = 0; i < dimM; i++)
		cout << " " << v[i];
	cout << endl;
}

void jacobi(int n, double a[][dimN], double b[], double x[]) {
	double sum, newX[dimN];
	int counter = 0;
	bool flag = false;
	cout << "Wyniki posrednie dla metody Jacobiego" << endl;
	cout << "i \t x[0] \t x[1]" << endl;
	while(1) {
		for (int i = 0; i < n; i++) {
			sum = 0;
			for (int j = 0; j < n; j++)
				if(j != i) sum += a[i][j] * x[j];
			newX[i] = (b[i] - sum) / a[i][i];
			if(fabs(x[i] - newX[i]) < eps) {
				flag = true;
				break;
			}
		}
		for (int i = 0; i < n; i++)
			x[i] = newX[i];
		cout << counter << "\t" << x[0] << "\t" << x[1] << "\t" << x[2] << endl;
		if(counter >= loop) {
			cout << "Osiagnieto ilosc iteracji." << endl;
			break;
		}
		if(flag) {
			cout << "Osiagnieto dokladnosc dla elementow wektora." << endl;
			break;
		}
		counter ++;
	}
}

void gaussSeidel(int n, double a[][dimN], double b[], double x[]) {
	double sum, previous;
	int counter = 0;
	bool flag = false;
	cout << "Wyniki posrednie dla metody Gaussa-Seidela" << endl;
	cout << "i \t x[0] \t x[1]" << endl;
	while(1) {
		for (int i = 0; i < n; i++) {
			sum = 0;
			for (int j = 0; j < n; j++)
				if(j != i) sum += a[i][j] * x[j];
			previous = x[i];
			x[i] = (b[i] - sum) / a[i][i];
			if(fabs(x[i] - previous) < eps) {
				flag = true;
				break;
			}
		}
		cout << counter << "\t" << x[0] << "\t" << x[1] << "\t" << x[2] <<endl;
		if(counter >= loop) {
			cout << "Osiagnieto ilosc iteracji." << endl;
			break;
		}
		if(flag) {
			cout << "Osiagnieto dokladnosc dla elementow wektora." << endl;
			break;
		}
		counter ++;
	}
}

void sor(int n, double a[][dimN], double b[], double x[], double omega) {
	double sum, newX;
	int counter = 0;
	bool flag = false;
	cout << "Wyniki posrednie dla metody SOR" << endl;
	cout << "i \t x[0] \t x[1]" << endl;
	while(1) {
		for (int i = 0; i < n; i++) {
			sum = 0;
			for (int j = 0; j < n; j++)
				if(j != i) sum += a[i][j] * x[j];
			newX = (1 - omega) * x[i] + omega *  (b[i] - sum) / a[i][i];
			if(fabs(x[i] - newX) < eps)
				flag = true;
			x[i] = newX;
			if(flag) break;
		}
		cout << counter << "\t" << x[0] << "\t" << x[1] << "\t" << x[2] <<endl;
		if(counter >= loop) {
			cout << "Osiagnieto ilosc iteracji." << endl;
			break;
		}
		if(flag) {
			cout << "Osiagnieto dokladnosc dla elementow wektora." << endl;
			break;
		}
		counter ++;
	}
}

int main(void) {
	double A[dimM][dimN], b[dimM], x[dimM], y[dimM];
	A[0][0] = 10.0; A[0][1] = 1.0; A[0][2] = 2.0;
	A[1][0] = 1.0; A[1][1] = 20.0; A[1][2] = 5.0;
	A[2][0] = 3.0; A[2][1] = 4.0; A[2][2] = 30.0;
	b[0] = 8.0; b[1] = -4.0; b[2] = -27.0;

	x[0] = 2.0; x[1] = 2.0; x[2] = 2.0;
	//Wypisanie macierzy
	cout << "Macierz wejsciowa" << endl;
	showMatrix(A);
	cout << endl;
	//Wypisanie wektora
	cout << "Wektor wejsciowy" << endl;
	showVector(b);
	cout << endl;
	//Wypisanie przyblizen poczatkowych
	cout << "Przyblizenia poczatkowe" << endl;
	showVector(x);
	cout << endl;
	//Wywolanie metod
	jacobi(dimM, A, b, x);
	cout << endl;
x[0] = 2.0; x[1] = 2.0; x[2] = 2.0;
	gaussSeidel(dimM, A, b, x);
	cout << endl;
x[0] = 2.0; x[1] = 2.0; x[2] = 2.0;
	sor(dimM, A, b, x, 0.5);
	cout << endl;
	system("pause");
	return 0;
}

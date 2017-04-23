#include <iostream>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <fstream>

using namespace std;

const double P = 1.0;
const double Q = 0.0;
const double R = 1.0;
const double S = 0.0;

const double ALPHA = 0.0;
const double BETA = 1.0;
const double GAMMA = -1.0;
const double FI = 0.0;
const double PSI = 1.0;
const double TETA = 0.0;

const double X1 = 0.0;
const double X2 = 2.0;

const int liczbaIteracji = 20;

//--------------------------------

double funkcjaAnalityczna(double x){
	return (exp(x / 2.0) - exp(2.0 - x / 2.0)) / (1.0 - exp(2.0));
}

void wektorB(double *b, int n){
    b[0] = -GAMMA;
    for(int i = 1; i < n - 1; i++)
        b[i] = S;
    b[n-1] = -TETA;
}

void metodaKonwencjonalna(double *wektorL, double *wektorD, double *wektorU, int n, double h){
    wektorU[0] = ALPHA / h;
    wektorD[0] = BETA - ALPHA / h;

    for(int i = 0; i < n - 2; i++)
	{
        wektorU[i+1] = P / (h * h) + Q / (2.0 * h);
        wektorD[i+1] = -2.0 * P / (h * h) + R;
        wektorL[i]   = P / (h * h) - Q / (2.0 * h);
    }

    wektorD[n-1] = FI / h + PSI;
    wektorL[n-2] = -FI / h;
}

void metodaNumerowa(double *wektorL, double *wektorD, double *wektorU, int n, double h){
    wektorU[0] = ALPHA / h;
    wektorD[0] = BETA - ALPHA / h;

    for(int i = 0; i < n - 2; i++){
        wektorU[i+1] = P / (h*h) + R / 12.0;
        wektorD[i+1] = -2.0 * P / (h * h) + R * 10.0 / 12.0;
        wektorL[i]   = P / (h * h) + R / 12.0;
    }

    wektorD[n-1] = FI / h + PSI;
    wektorL[n-2] = -FI / h;
}

void algThomasaMacierz(double *wektorL, double *wektorD, double *wektorU, int n){
    for(int i = 1; i < n; i++)
		wektorD[i] = wektorD[i] - wektorL[i-1] * wektorU[i-1] / wektorD[i-1];
}

void algThomasaWektor(double *wektorL,double *wektorD, double *wektorU, double *b, double *wynik, int n){
    for(int i = 1; i < n; i++)
		b[i] = b[i] - wektorL[i-1] * b[i-1] / wektorD[i-1];

	wynik[n-1] = b[n-1] / wektorD[n-1];

	for(int i = n - 2; i >= 0; i--)
		wynik[i] = (b[i] - (wektorU[i] * wynik[i+1])) / wektorD[i];
}

//--------------------------------

int main(){
	double h = 1.0, *wektorL, *wektorD, *wektorU, *b, *x, *wynik, xn, max, temp;
	int n;
	ofstream Dane;
	Dane.open("Dane.ods");

	cout << "log(h)\t m. Konwencjonalana\t M. Numerowa" << endl;

	for(int i = 0; i < liczbaIteracji; i++, h *= 0.5){
		n = (int) (X2 - X1) / h + 1;

        wektorL = new double[n-1];
        wektorD = new double[n];
        wektorU = new double[n-1];
        b = new double[n];
        wynik = new double[n];

        wektorB(b, n);
        metodaKonwencjonalna(wektorL, wektorD, wektorU, n, h);

        algThomasaMacierz(wektorL, wektorD, wektorU, n);
        algThomasaWektor(wektorL, wektorD, wektorU, b, wynik, n);

        max = 0.0;
		xn = X1;
		for(int i = 0; i < n; i++, xn += h){
			if(( temp = fabs(wynik[i] - funkcjaAnalityczna(xn))) > max)
				max = temp;
		}
		printf("%12g %12g", log10(h), log10(max));
		Dane << log10(h) << "\t";
		Dane << log10(max) << "\t";

		//-------------------------------

		wektorB(b, n);
		metodaNumerowa(wektorL, wektorD, wektorU, n, h);
        algThomasaMacierz(wektorL, wektorD, wektorU, n);
        algThomasaWektor(wektorL, wektorD, wektorU, b, wynik, n);

        max = 0.0;
		xn = X1;
		for(int i = 0; i < n; i++, xn += h){
			if(( temp = fabs(wynik[i] - funkcjaAnalityczna(xn))) > max)
				max = temp;
		}

		printf("%12g\n", log10(max));
		Dane << log10(max) << endl;

        delete []wektorL;
		delete []wektorD;
		delete []wektorU;
		delete []b;
		delete []wynik;
	}

	Dane.close();
	return 0;
}

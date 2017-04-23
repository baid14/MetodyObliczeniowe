#include "stdafx.h"
#include <sys/timeb.h>

void wykres_porownujacy(Projekt_D proj, double czas, int czestotliwosc, int metoda);
void wykres_max_blad(Projekt_D proj, int freq);
double Czas();
int ilosc = 400;
double czasStart;
double czasKoniec;


int main(){	
	Projekt_D proj;
	proj.informacje();
	proj.ustawZmienne(ilosc);
	proj.rozwAnalityczne();
	czasStart = Czas();
	//proj.rozw_laasonen(THOMAS);
	proj.rozwLaasonen(LU);

	//wykres_porownujacy(proj, proj.getDt(), 20, 1);
	///wykres_porownujacy(proj, proj.getDt(), 20, 2);
	
	//proj.computeKMB(ilosc);
	czasKoniec = Czas();
	//wykres_max_blad(proj, proj.getN()-1);
}


void wykresMaxBlad(Projekt_D proj, int freq)
{
		double dt = proj.getDt();
		int N = proj.getN();
		int M = proj.getM();

		double *pozycje;
		pozycje = proj.getPoz();

		double *czasy;
		czasy = proj.getCzas();

		fstream wyniki;
		wyniki.open("laasonen_KMB_wykres600.csv", ios::out);

		if(wyniki.is_open()){
			double **x;
			x = proj.getBezpo();
			
			double *ptrX;

			double **an;
			an = proj.getAnal();

			double err;
			double max;

			max = 0;
			ptrX = x[proj.getN()-1];
			for(int i = 0; i < M; i++){
				err = ptrX[i] - an[proj.getN()-1][i];
				err = fabs( err );
				if ( err > max ){
					max = err;
				}
			} 
				wyniki  << max << ";" << czasKoniec - czasStart << ";" << (proj.getA()/(double) (M-1)); 

			wyniki.close();
		}
		cout << "koniec" << endl;
		getchar();
}

double Czas(){
       struct timeb czas;
       double sekundy;
       ftime(&czas);

       sekundy= (double) czas.time;
       sekundy += (double) czas.millitm / 1000.0;
       return sekundy;
}
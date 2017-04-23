#include "stdafx.h"

void wykres_porownujacy(Projekt_C proj, double czas, int czestotliwosc, int metoda);
void wykres_max_blad(Projekt_C proj, int freq);

int main()
{
	Projekt_C proj;
	proj.informacje();

	proj.ustaw_zmienne(600);
	proj.rozw_anal();
	proj.rozw_laasonen(THOMAS);
	//proj.rozw_laasonen(LU);

	//wykres_porownujacy(proj, proj.getDt(), 20, 1);
	//wykres_porownujacy(proj, proj.getDt(), 20, 2);

	wykres_max_blad(proj, 50);
}

void wykres_porownujacy(Projekt_C proj, double czas, int czestotliwosc, int metoda)
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

}
#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

/*

Metody Obliczeniowe
By Wahuu

*/



double roz_analityczne(double t)
{
       return (1+t)*exp(-t);
}

// Metora Bezposrednia Euletra
double MBE(double t,double h,bool flaga)
{
       double i;
       double y = 1.0;
       double wart_dokl = 0.0;
       double blad = 0.0;
       for(i =0.0;i<t;i+=h)
       {
             wart_dokl = roz_analityczne(i);
             y = y + h * (-y+exp(-i));
             wart_dokl = fabs(wart_dokl - y);
		     if(wart_dokl > blad)
             blad = wart_dokl;
       }
       if(flaga == true)
                return y;
       else
           return blad;
}

// Metoda Posrednia Eulera
double MPE(double t,double h,bool flaga)
{
       double i;
       double y = 1.0;
       double wart_dokl = 0.0;
       double blad = 0.0;
       for(i =0.0;i<t;i+=h)
       {
             wart_dokl = roz_analityczne(i);
      	     y = (y + exp(-i) * h) / (1 + h);
      	     wart_dokl = fabs(wart_dokl - y);
		     if(wart_dokl > blad)
             blad = wart_dokl;
       }
        if(flaga == true)
                return y;
       else
           return blad;
}

// Metoda Trapezow
double MT(double t,double h,bool flaga)
{

       double i;
       double y = 1.0;
       double wart_dokl = 0.0;
       double blad = 0.0;
       for(i =0.0;i<t;i+=h)
       {
             wart_dokl = roz_analityczne(i);
      	     y = (2.0 * y -((y - exp(-i) - exp(-i + h)) * h)) / (2.0 + h);
      	     wart_dokl = fabs(wart_dokl - y);
		     if(wart_dokl > blad)
             blad = wart_dokl;
       }

       if(flaga == true)
                return y;
       else
           return blad;
}

int main()
{
    double krok = 0.1;
    double tmax = 2;
    fstream plik,bledy;

	plik.open("rozwiazania.csv", ios::out);
	if( plik.good() != true )
		cout << "Nie uzyskano dostepu do pliku!" << endl;
	plik << "t;Analityczne;MBE(stabline);MPE;MT"<<endl;
	for(double t = 0.0;t<=20;t+=0.01)
	{
         plik<<t<<";"<<roz_analityczne(t)<<";"<<MBE(t,0.01,true)<<";"<<MPE(t,0.01,true)<<";"<<MT(t,0.01,true)<<endl;
    }
    plik<<endl;
    for(double t = 0.0;t<=40;t+=4.0)
    {
               plik<<t<<";"<<roz_analityczne(t)<<";"<<MBE(t,4.0,true)<<endl;
    }
	bledy.open("bledy.csv", ios::out);
	if( bledy.good() != true )
		cout << "Nie uzyskano dostepu do pliku!" << endl;
	bledy << "krok;BME;PME;MT" << endl;
	for(int i=0; i<20; i++)
	{
		bledy <<  krok << ";" << MBE(tmax,krok,false)
			<<  ";" << MPE(tmax,krok,false) << ";" << MT(tmax ,krok,false) << endl;
		krok/=2;
	}
    //system("pause");
}

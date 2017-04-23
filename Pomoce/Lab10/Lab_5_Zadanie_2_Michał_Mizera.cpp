#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>

using namespace std;
bool zapis = false;
ofstream plik("dane10.txt");

double analityczne(double t)
{
	return 1 - exp( -10.0 * ( t + atan( t ) ));
}

double mbe(double n, double h)
{
	double y_nastepne = 0.0;
	double y  = 0.0;
	double bladMax = 0.0;
	double blad = 0.0;
	double t = h;

	for(double i = 0; i < n; i=i+1)
	{
		y_nastepne = y + h *( - ((10.0 *t * t + 20.0) /( t * t + 1.0) )*(y - 1.0) );
		double blad = fabs(y_nastepne - analityczne(t));

        if (zapis) plik<<t<<" "<<analityczne(t)<<" "<<y_nastepne<<endl;

		if (blad > bladMax) bladMax = blad;
		y = y_nastepne;
		t = t + h;
	}
	return bladMax;
}

double mpe(double n,double h)
{
	double y_nastepne = 0.0;
	double y  = 0.0;
	double bladMax = 0.0;
	double blad;
	double czynnik = 0.0;
    double t = h;

	for(double i = 0.0; i < n; i=i+1)
	{
	    czynnik = (10 * (t + h)*(t + h) + 20)/((t + h)*(t + h) + 1);
		y_nastepne = (y + h * czynnik )/( 1.0 + h* czynnik);
		double blad = fabs(y_nastepne - analityczne(t));

		if (zapis) plik<<t<<" "<<analityczne(t)<<" "<<y_nastepne<<endl;

		if (blad > bladMax) bladMax = blad;
		y = y_nastepne;
		t = t + h;
	}
	return bladMax;
}

double mt(double n,double h)
{
	double y_nastepne = 0.0;
	double y  = 0.0;
	double bladMax = 0.0;
	double blad;
	double czynnik1 = 0.0;
	double czynnik2 = 0.0;
    double t = h;

	for(double i = 0.0; i < n; i=i++)
	{
	    czynnik1 = (10.0 * t * t + 20.0)/(t * t + 1.0);
	    czynnik2 = (10.0 * (t + h)*(t + h) + 20.0)/((t + h)*(t + h) + 1.0);
		y_nastepne = (y - 0.5 * h * czynnik1 * (y - 1.0) + 0.5 * czynnik2 * h)/(1.0 + 0.5 * czynnik2 * h);
		double blad = fabs(y_nastepne - analityczne(t));

		if (zapis) plik<<t<<" "<<analityczne(t)<<" "<<y_nastepne<<endl;

		if (blad > bladMax) bladMax = blad;
		y = y_nastepne;
		t = t+h;
	}
	return bladMax;
}

int main()
{
    plik.precision(14);
	plik.flags(ios::fixed);

   // mbe(1000.0,0.01);
   // mpe(1000.0,0.01);
   // mt(1000.0,0.01);

    for ( double h = 1e-20; h < 0.0001 ; h*=2 )
     {
      if(!zapis)
      {
          plik<<log10(h)<<" "<<log10(mbe(1000.0,h))<<" "<<log10(mpe(1000.0,h))<<" "<<log10(mt(1000.0,h))<<endl;
      }
	}

    cout<<"Zakonczono obliczenia\n";
    return 0;
}

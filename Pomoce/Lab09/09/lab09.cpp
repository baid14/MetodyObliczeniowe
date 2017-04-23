/*#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <limits>
#include <iomanip>
using namespace std;

void wypiszWektor(double * wektorKroku, double *wektorBledu, int rozmiar)
{
	cout << fixed << setprecision(20);
	for (int i = 0; i < rozmiar; ++i) cout << wektorKroku[i] << "\t" << wektorBledu[i] << endl;
	cout << endl;
	std::cout.unsetf (std::ios_base::showbase);
	return;
}

double liczEpsilon ()
{
    double epsilon = 1.0, temp = epsilon + 1.0;
    while (temp > 1.0)
    {
		epsilon /= 2.0;
		temp = epsilon + 1.0;
	}
	return epsilon;
}

int liczIteracje ()
{
	int iteracje = 1;
	double epsilon = liczEpsilon();
	double h = 0.1;
	while (h > epsilon)
	{
		++iteracje;
		h /= 10.0;
	}
	return iteracje;
}

inline double metodaAnalityczna(double x)
{
	return (exp (x / 2.0) - exp (2.0 - x / 2.0)) / (1 - exp (2.0));
}

inline double *  rozwAnalityczne(double a, double b, double h, int rozmiar)
{
	double * x = new double [rozmiar],
	       * analityczne = new double [rozmiar];
	
	x [0] = a;
	x [rozmiar- 1] = b;
	
	analityczne[0] = metodaAnalityczna(x [0]);
	for (int i = 1; i < rozmiar - 1; ++i)
	{
		x [i] = x [i - 1] + h;
		analityczne[i] = metodaAnalityczna(x [i]);
	}
	analityczne[rozmiar - 1] = metodaAnalityczna(x [rozmiar - 1]);
	
	delete [] x;
	return analityczne;
}


inline void algorytmThomasa(double * l, double * d, double * u, double * x, double * b, int size)
{
	double * beta = new double [size];
	beta [0] = - u [0] / d [0];
	for (int i = 1; i < size; ++i)
		beta [i] = - u [i] / (l [i] * beta [i - 1] + d [i]);
	
	double * gamma = new double [size];
	gamma [0] = b [0] / d [0];
	for (int i = 1; i < size; ++i)
		gamma [i] = (b [i] - l [i] * gamma [i - 1]) / (l [i] * beta [i - 1] + d [i]);
	
	x [size - 1] = (b [size - 1] - l [size - 1] * gamma [size - 2]) / (l [size - 1] * beta [size - 2] + d [size - 1]);
	for (int i = size - 2; i >= 0; --i)
		x [i] = beta [i] * x [i + 1] + gamma [i];
	
	delete [] gamma;
	delete [] beta;
	return;
}

inline double maxBladKroku(double * numeryczne, double a, double b, double h, int rozmiar)
{
	double maxBlad = 0.0,
	       * analityczne = rozwAnalityczne(a, b, h, rozmiar);
	for (int i = 0; i <rozmiar; ++i)
		maxBlad = std::max (maxBlad, fabs (analityczne[i] - numeryczne[i]));
	delete [] analityczne;
	delete [] numeryczne;
	return maxBlad;
}

double liczTangens(double x_1, double x_2, double fx_1, double fx_2)
{
	return tan ((x_1 - x_2) / (fx_1 - fx_2));
}

// trzypunktowa dyskretyzacja konwencjonalna || metoda roznicowa
inline double * dyskretyzacjaKonw3p(double A, double B, double C, double h, double fa, double fb, int rozmiar)
{
	double * l = new double [rozmiar],
	       * d = new double [rozmiar],
	       * u = new double [rozmiar];
	
	u [rozmiar- 1] = l [0] = 0.0;
	d [rozmiar- 1] = d [0] = 1.0;
	
	double side_diagonal = A / (h * h), mid_diagonal = B - (2 * A / (h * h));
	for (int i = 1; i < rozmiar- 1; ++i) u [i - 1] = l [i] = side_diagonal;
	for (int i = 1; i < rozmiar- 1; ++i) d [i] = mid_diagonal;
	
	double * b = new double [rozmiar];
	b [0] = fa;
	b [rozmiar- 1] = fb;
	for (int i = 1; i < rozmiar- 1; ++i) b [i] = C;
	
	double * x = new double [rozmiar];
	
	algorytmThomasa(l, d, u, x, b, rozmiar);
	
	delete [] b;
	delete [] u;
	delete [] d;
	delete [] l;
	
	return x;
}

inline double * metodaNumerowa(double A, double B, double C, double h, double fa, double fb, int rozmiar)
{
	double * l = new double [rozmiar];
	double * d = new double [rozmiar];
	double * u = new double [rozmiar];
	u [rozmiar - 1] = l [0] = 0.0;
	d [rozmiar - 1] = d [0] = 1.0;
	
	double mid_diagonal = (24.0 * A - 10.0 * B * h * h) / (12.0 * A - B * h * h);
	for (int i = 1; i < rozmiar - 1; ++i) u [i - 1] = l [i] = -1.0;
	for (int i = 1; i < rozmiar - 1; ++i) d [i] = mid_diagonal;
	
	double * b = new double [rozmiar];
	b [0] = fa;
	b [rozmiar- 1] =fb;
	for (int i = 1; i < rozmiar- 1; ++i)
		b [i] = C;
	double * x = new double [rozmiar];
	
	algorytmThomasa(l, d, u, x, b, rozmiar);
	delete [] b;
	delete [] u;
	delete [] d;
	delete [] l;
	
	return x;
}
void zapiszDoPliku(double *wektorKrokow, double *wektorBledu, int ile_krokow, string nazwa);

void metodaNumeryczna(double *(*funkcja) (double , double , double , double , double , double , int ), string nazwa)
{
	int ileKrokow = liczIteracje();                                // ileKrokow = 16
	int rozmiar;
	double a = 0.0, fa = 1.0,
	       b = 2.0, fb = 0.0,
	       h = 0.1,
	       A = 1.0, B = -0.25, C = 0.0;
	double * wektorBledu = new double [ileKrokow],
	       * wektorKrokow = new double [ileKrokow];
	
	for (int i = 0; i < ileKrokow; ++i)
	{
		rozmiar = static_cast <int> ((fabs (b - a) / h) + 1.0);
		wektorBledu [i] = log10 (maxBladKroku(funkcja(A, B, C, h, fa, fb, rozmiar), a, b, h, rozmiar));
		wektorKrokow[i] = log10 (h);
		h /= 2.0;
	}
	

	//for(int i=0; i<ileKrokow-1; i++)
	//cout << "Rzad dokladnosci wynosi: " << liczTangens(wektorKrokow[i], wektorKrokow[i+1], wektorBledu[i], wektorBledu[i+1]) << endl;
	cout << "Rzad dokladnosci wynosi: " << liczTangens(wektorKrokow[0], wektorKrokow[1], wektorBledu[0], wektorBledu[1]) << endl;
	//wypiszWektor(wektorKrokow, wektorBledu, ileKrokow);
	zapiszDoPliku(wektorKrokow, wektorBledu, ileKrokow, nazwa);

	delete [] wektorBledu;
	return;
}

int main (void)
{
    for(int i=0; i<10; i++) {
	metodaNumeryczna(dyskretyzacjaKonw3p, "metoda-dyskretyzacja");}
	//metodaNumeryczna(metodaNumerowa, "metoda-numerowa");
	
	system("PAUSE");
	return 0;
}

void zapiszDoPliku(double *wektorKrokow, double *wektorBledu, int ile_krokow, string nazwa)
{
	fstream plik;
	string tytul=nazwa;
	tytul+=".csv";
	plik.open(tytul.c_str(),fstream::out);
	for (int i=0; i<ile_krokow; i++) plik << wektorKrokow[i] <<","<< wektorBledu[i] << endl;
	cout << "Zapisano!\n\n";
	plik.close();
}*/



#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

double liczEpsilon ()
{
    double epsilon = 1.0, temp = epsilon + 1.0;
    while (temp > 1.0)
    {
		epsilon /= 2.0;
		temp = epsilon + 1.0;
	}
	return epsilon;
}

int liczIteracje ()
{
	int iteracje = 1;
	double epsilon = liczEpsilon();
	double h = 0.1;
	while (h > epsilon)
	{
		++iteracje;
		h /= 10.0;
	}
	return iteracje;
}

void algThomasa(double *l, double *d, double *u, int rozmiar)
{	                                                          
    	 for(int i=1; i<rozmiar; i++) d[i]=d[i]-((l[i-1]/d[i-1])*u[i-1]);
}

void rozwThomas(double *l, double *d, double *u, double *b, double *x, int rozmiar)
{
	for(int i=1; i<rozmiar; i++) b[i]=b[i]-((l[i-1]/d[i-1])*b[i-1]);  

	//Zabezpieczenie przed dzieleniem przez 0
    	if(d[rozmiar-1]==0) b[rozmiar-1]=0;
    	else x[rozmiar-1]=b[rozmiar-1]/d[rozmiar-1];

    	for(int i=rozmiar-2; i>=0; i--) x[i]=(b[i]-u[i]*x[i+1])/d[i];
}

void analityczne(double *wektor, double h, int n) 
{    
	 //obliczanie analitycznych rozwiazan
     double krok=h;
	 //Pierwszy el wektora wynikowego jest równy wartoœci funkcji w punkcie startowym
     wektor[0]=1.0;
     for(int i=1; i<=n-1; i++) 
	 {
		 wektor[i]=(exp(h)-exp(2.0-h))/(1.0-exp(2.0));
         h+=krok;
     }
	 //Ostatni el wektora wynikowego jest równy wartoœci funkcji w punkcie koñcowym
	 wektor[n]=0.0;
}
double *metodaKonwencjonalna(double h, int n) 
{ 
	double *l=new double[n]; //deklaracja wektorow do przechowywania
	double *d=new double[n+1]; //przekatnych macierzy trojdiagonalnej
	double *u=new double[n]; //wektorow rozwiazani i bledow
	double *b=new double[n+1];
	double *x=new double[n+1];
	//Wype³niamy wektory l, d, u, b
	for(int i=0; i<=n; i++) 
	{
     	l[i]=u[i]=-1.0;
        d[i]=(h*h + 2.0);
        b[i]=0.0;
    }
	//usuniêty d[n]
	d[0]=b[0]=1.0;
	u[0]=l[n-1]=0.0;

    algThomasa(l,d,u,n);
    rozwThomas(l,d,u,b,x,n);

	delete [] b;
	delete [] u;
	delete [] d;
	delete [] l;
	return x;
}

double *metodaNumerowa(double h, int n) 
{ 
	double *l=new double[n]; //deklaracja wektorow do przechowywania
	double *d=new double[n+1]; //przekatnych macierzy trojdiagonalnej
	double *u=new double[n]; //wektorow rozwiazani i bledow
	double *b=new double[n+1];
	double *x=new double[n+1];
	//wype³niamy wektory l, d, u, b
    for(int i=0; i<=n; i++) 
	{
		l[i]=u[i]=h*h-12.0;
        d[i]=24.0+10.0*h*h;
        b[i]=0.0;
    }
    //usuniêty d[n]
	d[0]=b[0]=1.0;
	u[0]=l[n-1]=0.0;
    algThomasa(l,d,u,n);
    rozwThomas(l,d,u,b,x,n);
	delete [] b;
	delete [] u;
	delete [] d;
	delete [] l;
	return x;
}

double maxBlad(double *rozw, double *b, double *blad, int n)
{
	//liczymy wektor b³êdu metody
	for(int i=0; i<n-1; i++) blad[i]=b[i+1]-rozw[i+1];
	double max=0.0;
	//szukamy max b³êdu w wektorze
	for(int i=0; i<n-1; i++) if(max<fabs(blad[i])) max=fabs(blad[i]);
    //cout<<max<<endl;
	return max;
}

int main() 
{
	double h,g;
    int e,nr;
    cout<<"Ktory schemat chcesz uzyc?"<<endl;
    cout<<"1. Schemat konwencjonalny"<<endl;
    cout<<"2. Schemat Numerowa"<<endl;
    cin>>nr;
    cout<<"podaj krok h"<<endl;
    cin>>h;
	
	int ile=liczIteracje();
	cout << "ile=" << ile << endl;
	
    switch(nr) {
      	case 1: 
		{
			double *blad2;
			double *rozwAnalityczne;
			double *rozwNumeryczne;
			for(int i=0; i<ile; i++) 
			{
			g=(1.0/h); //liczymy liczbe przedzialow
			e=(int)g;
			rozwAnalityczne=new double[e+1];
			rozwNumeryczne=new double[e+1];
			blad2=new double[e+1];
			analityczne(rozwAnalityczne,h,e);
            rozwNumeryczne=metodaKonwencjonalna(h,e);
            cout<<"norma maksimum z wektora bledow wynosi: " << maxBlad(rozwAnalityczne,rozwNumeryczne,blad2,e) <<endl;
			h*=0.5;
			}
			delete rozwNumeryczne;
			delete rozwAnalityczne;
			delete blad2;

        } break;
        case 2: 
		{
			for(int i=0; i<ile; i++) {
			g=(1.0/h); //liczymy liczbe przedzialow
			e=(int)g;
			double *rozwAnalityczne=new double[e+1];
			double *rozwNumeryczne=new double[e+1];
			double *blad2=new double[e+1];
			analityczne(rozwAnalityczne,h,e);
            rozwNumeryczne=metodaNumerowa(h,e);
            cout<<"norma maksimum z wektora bledow wynosi:"<< maxBlad(rozwAnalityczne,rozwNumeryczne,blad2,e) <<endl;
			h*=0.5;}
        }
	}
    system("PAUSE");
    return 0;
}

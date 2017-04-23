#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>

using namespace std;

ofstream plik("dane.txt");

double analityczne(double x)
{
	return (exp(x / 2.0) - exp(2.0 - x / 2.0)) / (1.0 - exp(2.0));
}
double norma_maksimum(double *w, int n)
{
	double maximum = fabs(w[0]);
	for(int i=1; i<n; i++)
		if(fabs(w[i] > maximum))
			maximum = fabs(w[i]);
	return maximum;
}

void thomas_elim(double *a, double *d, double *c, double *f, int n)
{
    for(int i = 1;i < n; i++)
    {
        f[i-1] = (a[i-1] / d[i-1]); // wspolczynnik
        d[i]=d[i] - f[i-1] * c[i -1];
        a[i-1] = 0;
    }
}

void thomas_rown(double *d,double *c,double *b,double* f,double *x, int n)
{
    for(int i = 1;i < n ; i++)
      b[i] = b[i] - f[i-1]*b[i-1];

    x[n-1] = b[n-1]/d[n-1];

    for(int i = n - 2; i >= 0 ; i--)
    {
        x[i] = (b[i] - c[i] * x[i+1])/d[i];
    }
}

template<typename typ>
void wyswietlWektor(typ *wektor,int n)
{
    cout<<endl;
	for (int i = 0; i<n; i++)
      printf("%9.2lf",wektor[i]);
	cout<<endl;
}

template <typename typ>
typ* alokujWektor(int n)
{
    typ *wektor;
    wektor = new typ[n];
    return wektor;
}

void usunWektor(double *wektor)
{
    delete []wektor;
}

double dyskretyzacja_konwencjonalna(int n, double xpocz, double xkon,double* wektor_wynikowy)
{
	double punkt=xpocz;
	double *a = alokujWektor<double>(n - 1); // ponizej diagonalnej
	double *b = alokujWektor<double>(n); // wektor b
	double *c = alokujWektor<double>(n - 1); // powyzej diagonalnej
	double *d = alokujWektor<double>(n); // diagonalna
	double *blad = alokujWektor<double>(n); // wektor bledu
	double *x = alokujWektor<double>(n);	 // rozwiazania
    double *f = alokujWektor<double>(n); // wektor wspolczynnikow z eliminacji
    double h = fabs(xkon - xpocz)/(n-1);

    d[0] = 1.0;
	b[0] = 1.0;
	c[0] = 0.0;

	for(int i=1; i<n-1; i++)
    {
		d[i] = -(8.0+h*h)/(4.0*h*h);
		b[i] = 0.0;
	}

    for(int i = 0; i < n-2 ; i++)
    {
        a[i] = 1.0/(h*h);
    }

    for(int i = 1; i < n-1 ; i++)
    {
        c[i] = 1.0/(h*h);
    }

	b[n-1] = 0.0;
	d[n-1] = 1.0;
	a[n-2] = 0.0;

    thomas_elim(a,d,c,f,n);
    thomas_rown(d,c,b,f,x,n);

//  cout<<setw(15)<<"liczba iteracji:"<<setw(20)<<"krok h:"<<setw(25)<<"wartosc przyblizona"<<setw(20)<<"wartosc dokladna"<<endl;

	for(int i=0; i<n; i++)
    {
//      cout<<setw(15)<<i<<setw(20)<<punkt<<setw(25)<<x[i]<<setw(20)<<analityczne(punkt)<<endl;
//      plik<<punkt<<" "<<analityczne(punkt)<<" "<<x[i]<<" "<<endl;
         wektor_wynikowy[i] = x[i];
         blad[i] = fabs(x[i] - analityczne(punkt));
		 punkt += h;
	}
    usunWektor(a),usunWektor(b),usunWektor(c),usunWektor(d),usunWektor(x),usunWektor(f);

	return norma_maksimum(blad, n);
}

double dyskretyzacja_numerowa(int n, double xpocz, double xkon,double* wektor_wynikowy)
{
	double punkt=xpocz;
	double *a = alokujWektor<double>(n-1); // ponizej diagonalnej
	double *b = alokujWektor<double>(n); // wektor b
	double *c = alokujWektor<double>(n-1); // powyzej diagonalnej
	double *d = alokujWektor<double>(n); // diagonalna
	double *blad = alokujWektor<double>(n); // wektor bledu
	double *x = alokujWektor<double>(n);	 // rozwiazania
    double *f = alokujWektor<double>(n); // wektor wspolczynnikow z eliminacji
    double h = fabs(xkon - xpocz)/(n-1);

    d[0] = 1.0;
	b[0] = 1.0;
	c[0] = 0.0;

    for(int i=1; i<n-1; i++)
    {
		d[i] = -(48.0 + 5.0*h*h)/(24.0*h*h);
		b[i] = 0.0;
	}

    for(int i = 0; i < n-2 ; i++)
    {
       a[i] = (48.0 - h*h)/(48.0*h*h);
    }

    for(int i = 1; i < n-1 ; i++)
    {
       c[i] = (48.0 - h*h)/(48.0*h*h);
    }

	b[n-1] = 0.0;
	d[n-1] = 1.0;
	a[n-2] = 0.0;

    thomas_elim(a,d,c,f,n);
    thomas_rown(d,c,b,f,x,n);

//  cout<<setw(15)<<"liczba iteracji:"<<setw(20)<<"krok h:"<<setw(25)<<"wartosc przyblizona(numerow)"<<setw(20)<<"wartosc dokladna"<<endl;

	for(int i=0; i<n; i++)
    {
//       cout<<setw(15)<<i<<setw(20)<<punkt<<setw(25)<<x[i]<<setw(20)<<analityczne(punkt)<<endl;
//       plik<<punkt<<" "<<analityczne(punkt)<<" "<<x[i]<<" "<<endl;
         wektor_wynikowy[i] = x[i];
         blad[i] = fabs(x[i] - analityczne(punkt));
		 punkt += h;
	}

    usunWektor(a),usunWektor(b),usunWektor(c),usunWektor(d),usunWektor(x),usunWektor(f);

	return norma_maksimum(blad, n);
}

void obliczaj(int n, int nmax, double xpocz, double xkon)
{
	double h;
    double* konwencjonalna_wyniki = alokujWektor<double>(nmax);
    double* numerowa_wyniki = alokujWektor<double>(nmax);
//	cout<<setw(10)<<"N: "<<setw(20)<<"krok:"<<setw(20)<<"blad konwencjonalnej: "<<setw(20)<<"blad numerowa"<<endl;
	while(n <= nmax)
    {
        h = fabs(xkon - xpocz)/(n-1);
//      cout<<setw(10)<<n<<setw(20)<<h<<setw(20)<<dyskretyzacja_konwencjonalna(n,xpocz,xkon,)<<endl;
        plik<<log10(h)<<" "<<log10(dyskretyzacja_konwencjonalna(n,xpocz,xkon,konwencjonalna_wyniki))<<" "<<log10(dyskretyzacja_numerowa(n,xpocz,xkon,numerowa_wyniki))<<endl;
        n = n+10;
	}
	usunWektor(konwencjonalna_wyniki),usunWektor(numerowa_wyniki);
}

void rysuj_do_pliku(int n,double xpocz,double xkon)
{
  double* konwencjonalna_wyniki = alokujWektor<double>(n);
  double* numerowa_wyniki = alokujWektor<double>(n);
  double punkt = xpocz;

  dyskretyzacja_konwencjonalna(n,xpocz,xkon,konwencjonalna_wyniki);
  dyskretyzacja_numerowa(n,xpocz,xkon,numerowa_wyniki);

   double h = fabs(xkon - xpocz)/(n-1);

  for(int i = 0 ; i < n; i++)
  {
      plik<<punkt<<" "<<konwencjonalna_wyniki[i]<<" "<<numerowa_wyniki[i]<<" "<<analityczne(punkt)<<endl;
      punkt += h;
  }
  usunWektor(konwencjonalna_wyniki),usunWektor(numerowa_wyniki);
}

int main()
{
    plik.precision(14);
	plik.flags(ios::fixed);

     obliczaj(10,2000,0.0,2.0);
     //rysuj_do_pliku(50,0,2.0);

    cout<<"Zakonczono obliczenia\n";
    return 0;
}

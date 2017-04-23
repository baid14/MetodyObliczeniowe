#include<iostream>
#include<cmath>
using namespace std;

double fun(double x)
{ 
	return exp(x);
}
double centralnasr(double a, double h)              //f-cja obliczajaca roznice centralna dla sr wartosci
{
	double f;
	f = (fun(a+h) - fun(a-h))/(2*h);
	return f;	
}

double wstecznasr(double a, double h)			//f-cja obliczajaca roznice wsteczna dla sr wartosci
{
	double f;
	f = (fun(a)-fun(a-h))/(h);
	return f;
}

double progresywnasr (double a, double h)		//f-cja obliczajaca roznice progresywna dla wart srodkowej
{
	double f;
	f = (fun(a+h) - fun(a))/h;
	return f;
}

double kr01(double a, double h)				//f-cja obliczajaca rozniczke dla wartosci poczatkowej
{
	double f;
	f = (-3./2.*fun(a) + 2*fun(a+h) - 1./2.*fun(a+2*h))/h;
	return f;

}

double kr02(double a, double h)				//f-cja obliczajaca rozniczke dla wartosci koncowej
{
	double f;
	f = (1.0/2.0*fun(a-2*h) - 2.0*fun(a-h) + 3.0/2.0*fun(a))/h;
	return f;

}

double blad(double a, double b)				//f-cja obliczajaca blad			
{
	double c;
	c = fabs((a-b)/b);
	return c;
}


int main()
{
	double a=0., b=1.;
	double c;
	double k,l,m,n,t,r,w,q,y,z;
	c = (a+b)/2;
	double h=0.1;
	double st;  		//zmienna okreslajaca warunek stopu
cout << "roznica centralna: " << centralnasr(c,h) << "\t" << "roznica wsteczna: " << wstecznasr(c,h) << "\t" << "roznica progresywna: " << progresywnasr(b,h) << endl;
cout << "wartosc rzeczywista pochodnej: " << exp(c) << endl;
	cout << "dla wartosci srodkowej najlepsza pochodna liczona jest z roznicy centralnej,wiec wykorzystuje ja w pozniejszych oblicezniach" << endl;
	cout << "logarytm bledu pocz: " << " | " << "logarytm bledu srodka: " << " | " << "logarytm bledu koncowego: " << " | " << "logarytm h" << endl;
	while(1)
	{
	k = kr01(a,h);               //zmienna przechowujaca wart pochodnej dla poczatkowego kranca
	l = centralnasr(c,h);		//zmienna przech wart pochodnej dla punktu srodkowego
	m = kr02(b,h);			// zmienna przech wart pochodnej dla punktu koncowego
	n = fun(a);			
	t = fun(b);			//zmienne przechowujace rzeczywiste wart w zadanych pkt
	r = fun(c);
	
	cout << "blad dla wartosci poczatkowej: " << blad(k,n) << "\t";
	cout << "blad dla wart srodkowej: " << blad(l,r) << "\t";
	cout << "blad dla wart koncowej: " << blad(m,t) << endl;
	cout << "logarytm bledu pocz: " << " | " << "logarytm bledu srodka: " << " | " << "logarytm bledu koncowego: " << " | " << "logarytm h" << endl;
	w = log10(blad(k,n));                 
	q = log10(blad(l,r));
	y = log10(blad(m,t));
	z = log10(h);
	cout << "\t" <<  w << "\t|\t" << q << "\t|\t" << y << "\t\t|\t" << z << endl;
	st = h+1.0;
	if (st==1.0) break; 
	h/=10;
	cout << endl << endl;
	}
	system("PAUSE");
	return 0;  
	


}

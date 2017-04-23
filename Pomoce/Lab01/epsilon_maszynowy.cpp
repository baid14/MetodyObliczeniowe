#include <iostream>
using namespace std;


class Epsilon 
{
public:
	float licz_F(float baza);
	double licz_D(double baza);
	long double licz_L(long double baza);
	Epsilon() {}
	~Epsilon() {}
};

//double
double Epsilon::licz_D(double baza)
{
	int liczba_bman;
	double eps, dxpl;
	liczba_bman = 0; eps=0.5; dxpl = baza+eps;
	while (dxpl>baza)
	{
		eps = eps/2;
		dxpl = baza+eps;
		liczba_bman++;
	}
	cout << "Liczba bitow mantysys (double) " << liczba_bman << endl <<"Epsilon maszynowy: ";
	return eps;
}


//long double
long double Epsilon::licz_L(long double baza)
{
	int liczba_bman;
	long double eps, dxpl;
	liczba_bman=0; eps=0.5; dxpl = baza+eps;
	while (dxpl>baza)
	{
		eps = eps/2;
		dxpl = baza+eps;
		liczba_bman++;
	}
	cout << "Liczba bitow mantysys (long) " << liczba_bman << endl <<"Epsilon maszynowy: ";
	return eps;
}

//float
float Epsilon::licz_F(float baza)
{
	int liczba_bman;
	float eps, dxpl;
	liczba_bman=0; eps=0.5; dxpl = baza+eps;
	while (dxpl>baza)
	{
		eps = eps/2;
		dxpl = baza+eps;
		liczba_bman++;
	}
	cout << "Liczba bitow mantysys (float) " << liczba_bman << endl <<"Epsilon maszynowy: ";
	return eps;
}



int main()
{
	Epsilon eps;
	cout << "Baza = 1" << endl << endl << "float" << endl;
	cout << eps.licz_F(1.0) << endl << endl;
	cout << "double" << endl;
	cout << eps.licz_D(1.0) << endl << endl;
	cout << "long double" << endl;
	cout << eps.licz_L(1.0) << endl << endl;

	cout << "Baza = 100" << endl << endl << "float" << endl;
	cout << eps.licz_F(100.0) << endl << endl;
	cout << "double" << endl;
	cout << eps.licz_D(100.0) << endl << endl;
	cout << "long double" << endl;
	cout << eps.licz_L(100.0) << endl << endl;

    cout << "Baza = 1e4" << endl << endl << "float" << endl;
	cout << eps.licz_F(1e4) << endl << endl;
	cout << "double" << endl;
	cout << eps.licz_D(1e4) << endl << endl;
	cout << "long double" << endl;
	cout << eps.licz_L(1e4) << endl << endl;

    cout << "Baza = 1e8" << endl << endl << "float" << endl;
	cout << eps.licz_F(1e8) << endl << endl;
	cout << "double" << endl;
	cout << eps.licz_D(1e8) << endl << endl;
	cout << "long double" << endl;
	cout << eps.licz_L(1e8) << endl << endl;

return 0;
}


#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <typeinfo>

using namespace std;

int krokow = 15;
int metod = 7;

template <typename T> T roznicaProgresywna2pkt(T x_0, T x_plus);
template <typename T> T roznicaWsteczna2pkt(T x_minus, T x_0);
template <typename T> T roznicaCentralna2pkt(T x_minus, T x_plus);
template <typename T> T roznicaProgresywna3pkt(T x_0, T x_plus1, T x_plus2);
template <typename T> T roznicaWsteczna3pkt(T x_minus2, T x_minus1, T x_0);
template <typename T> T pochodna(T x);
template <typename T> void obliczRoznice(T * tablicaKrokow, T ** tablicaWynikow);
template <typename T> void zapiszDoPliku(string * tablicaOpisow, T * tablicaKrokow, T ** tablicaWynikow);

int main(int argc, char * argv[])
{	
	float * krokiFloat;
	double * krokiDouble;
	long double * krokiLongDouble;
	
	krokiFloat = new float [krokow];
	krokiDouble = new double [krokow];
	krokiLongDouble = new long double [krokow];
	
	float ** wynikiFloat;
	double ** wynikiDouble;
	long double ** wynikiLongDouble;
	
	wynikiFloat = new float * [metod];
	wynikiDouble = new double * [metod];
	wynikiLongDouble = new long double * [metod];
	
	for(int i = 0; i < metod; i++)
	{
		wynikiFloat[i] = new float [krokow];
		wynikiDouble[i] = new double [krokow];
		wynikiLongDouble[i] = new long double [krokow];
	}
	
	string * opisy;
	
	opisy = new string [metod];
	
	opisy[0] = "Punkt środkowy - różnica progresywna dwupunktowa";
	opisy[1] = "Punkt środkowy - różnica wsteczna dwupunktowa";
	opisy[2] = "Punkt środkowy - różnica centralna dwupunktowa";
	opisy[3] = "Punkt środkowy - różnica progresywna trzypunktowa";
	opisy[4] = "Punkt środkowy - różnica wsteczna trzypunktowa";
	
	// opisy[5] = "Punkt końcowy - różnica progresywna dwupunktowa";
	opisy[5] = "Punkt końcowy - różnica wsteczna dwupunktowa";
	// opisy[7] = "Punkt końcowy - różnica centralna dwupunktowa";
	// opisy[8] = "Punkt końcowy - różnica progresywna trzypunktowa";
	opisy[6] = "Punkt końcowy - różnica wsteczna trzypunktowa";
	
	obliczRoznice(krokiFloat, wynikiFloat);
	obliczRoznice(krokiDouble, wynikiDouble);
	obliczRoznice(krokiLongDouble, wynikiLongDouble);
	
	zapiszDoPliku(opisy, krokiFloat, wynikiFloat);
	zapiszDoPliku(opisy, krokiDouble, wynikiDouble);
	zapiszDoPliku(opisy, krokiLongDouble, wynikiLongDouble);
	
	return 0;
}

// Dwupunktowe...
template <typename T> T roznicaProgresywna2pkt(T x_0, T x_plus)
{
	return (cos(x_plus) - cos(x_0)) / (x_plus - x_0);
}

template <typename T> T roznicaWsteczna2pkt(T x_minus, T x_0) 
{
	return (cos(x_0) - cos(x_minus)) / (x_0 - x_minus);
}

template <typename T> T roznicaCentralna2pkt(T x_minus, T x_plus) 
{
	return (cos(x_plus) - cos(x_minus)) / (x_plus - x_minus);
}

// Trzypunktowe...
template <typename T> T roznicaProgresywna3pkt(T x_0, T x_plus1, T x_plus2)
{
	return (-3 * cos(x_0) + 4 * cos(x_plus1) - cos(x_plus2)) / (x_plus2 - x_0);
}

template <typename T> T roznicaWsteczna3pkt(T x_minus2, T x_minus1, T x_0)
{
	return (3 * cos(x_0) - 4 * cos(x_minus1) + cos(x_minus2)) / (x_0 - x_minus2);
}

// Prawdziwa pochodna...
template <typename T> T pochodna(T x)
{
	return -sin(x); 
}

// Zapisujemy w tablicy wyniki dla różnych metod i różnych kroków.
template <typename T> void obliczRoznice(T * tablicaKrokow, T ** tablicaWynikow)
{
	T b = M_PIl / 4.0;
	T c = M_PIl / 2.0;
	
	T p;
	
	T h = 0.1;
	
	for (int i = 0; i < krokow; i++)
	{
		// Pokazujemy...
		cout << "Obliczam różnice (typ = " << typeid(T).name() << ", krok = " << h << ")..." << endl;
		
		// Zapisujemy krok:
		tablicaKrokow[i] = h;
		
		// Obliczamy:
		p = pochodna(b);
		
		tablicaWynikow[0][i] = fabs(p - roznicaProgresywna2pkt(b, b + h));
		tablicaWynikow[1][i] = fabs(p - roznicaWsteczna2pkt(b - h, b));
		tablicaWynikow[2][i] = fabs(p - roznicaCentralna2pkt(b - h, b + h));
		tablicaWynikow[3][i] = fabs(p - roznicaProgresywna3pkt(b, b + h, b + 2 * h));
		tablicaWynikow[4][i] = fabs(p - roznicaWsteczna3pkt(b - 2 * h, b - h, b));
		
		p = pochodna(c);
		
		// tablicaWynikow[5][i] = fabs(p - roznicaProgresywna2pkt(c, c + h));
		tablicaWynikow[5][i] = fabs(p - roznicaWsteczna2pkt(c - h, c));
		// tablicaWynikow[7][i] = fabs(p - roznicaCentralna2pkt(c - h, c + h));
		// tablicaWynikow[8][i] = fabs(p - roznicaProgresywna3pkt(c, c + h, c + 2 * h));
		tablicaWynikow[6][i] = fabs(p - roznicaWsteczna3pkt(c - 2 * h, c - h, c));
		
		// Zmniejszamy krok:
		h = h * 0.1;
	}
}

template <typename T> void zapiszDoPliku(string * tablicaOpisow, T * tablicaKrokow, T ** tablicaWynikow)
{
	fstream f;
	string nazwa = "Dane - typ ";
	
	nazwa += typeid(T).name();
	nazwa += ".csv";
	
	cout << "Zapisuje..." << endl;
	
	f.open(nazwa.c_str(), fstream::out);
	
	// f << "Typ: " << typeid(T).name() << endl << endl;
	
	f << "\"Krok\",";
	
	for(int i = 0; i < metod; i++)
	{
		f << "\"" << tablicaOpisow[i] << "\",";
	}
	
	f << endl;
	
	for(int i = 0; i < krokow; i++)
	{
		f << tablicaKrokow[i] << ",";
		
		for(int j = 0; j < metod; j++)
		{
			f << tablicaWynikow[j][i] << ",";
		}
		
		f << endl;
	}
	
	f.close();
	system("PAUSE");
}

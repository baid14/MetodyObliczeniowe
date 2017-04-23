// Program implementujacy algorytm Thomasa.

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>

using namespace std;

class Exception
{
	protected:
		string msg;
	public:
		Exception() { msg = ""; }
		Exception(string msg) { this->msg = msg; }
		string getMsg() { return msg; }
};

class ConversionException : public Exception
{
	public:
		ConversionException() : Exception() { };
		ConversionException(string msg) : Exception(msg) { };
};

class SolveException : public Exception
{
	public:
		SolveException() : Exception() { };
		SolveException(string msg) : Exception(msg) { };
};

template <class T> class Converter
{
	public:
		static T fromString(const string & s)
		{
			istringstream stream(s);
			T val;
			
			if((stream >> noskipws >> val).fail())
			{
				throw ConversionException("Podany ciag nie moze zostac przekonwertowany na liczbe!");
			}
			
			return val;
		}
};

void printVector(long double * elem, int dim)
{
	cout << setprecision(5) << fixed;
	
	cout << "[";
	for(int i = 0; i < dim; i++)
	{
		cout << setw(10) << elem[i] << " ";
	}
	cout << "]" << endl;
}

void printTridiagonalMatrix(long double * sup, long double * dia, long double * sub, int dim)
{
	cout << setprecision(5) << fixed;
	
	for(int i = 0; i < dim; i++)
	{
		cout << "|";
		if(i == 0)
		{
			cout << setw(10) << dia[i] << " ";
			cout << setw(10) << sup[i] << " ";
			for(int j = 2; j < dim; j++)
			{
				cout << "           ";
			}
		}
		else if (i == dim - 1)
		{
			for(int j = 0; j < dim - 2; j++)
			{
				cout << "           ";
			}
			cout << setw(10) << sub[i - 1] << " ";
			cout << setw(10) << dia[i] << " ";
		}
		else
		{
			for(int j = 0; j < i - 1; j++)
			{
				cout << "           ";
			}
			cout << setw(10) << sub[i - 1] << " ";
			cout << setw(10) << dia[i] << " ";
			cout << setw(10) << sup[i] << " ";
			for(int j = i + 2; j < dim; j++)
			{
				cout << "           ";
			}
		}
		
		cout << "|" << endl;
	}
}

int main(int argc, char* argv[])
{
	string dimString;
	string elString;
	int dim;
	long double * elemSub;
	long double * elemDia;
	long double * elemSup;
	long double * elemB;
	long double * x;
	
	/*
	 * Wczytujemy dane.
	 */
	
	dim = 5;
	
	elemSub = new long double [dim - 1];
	elemDia = new long double [dim];
	elemSup = new long double [dim - 1];
	elemB = new long double [dim];
	x = new long double [dim];
	
	elemSub[0] = 1.0;
	elemSub[1] = 2.0;
	elemSub[2] = 3.0;
	elemSub[3] = 4.0;
	
	elemDia[0] = 10.0;
	elemDia[1] = 11.0;
	elemDia[2] = 12.0;
	elemDia[3] = 13.0;
	elemDia[4] = 14.0;
	
	elemSup[0] = 4.0;
	elemSup[1] = 3.0;
	elemSup[2] = 2.0;
	elemSup[3] = 1.0;
	
	elemB[0] = 18.0;
	elemB[1] = 32.0;
	elemB[2] = 44.0;
	elemB[3] = 36.0;
	elemB[4] = 22.0;
	
	
	
	cout << "Podany wektor poddiagonalny:" << endl;
	
	printVector(elemSub, dim - 1);
	
	cout << "Podany wektor diagonalny:" << endl;
	
	printVector(elemDia, dim);
	
	cout << "Podany wektor naddiagonalny:" << endl;
	
	printVector(elemSup, dim - 1);
	
	cout << "Podany wektor wyrazow wolnych:" << endl;
	
	printVector(elemB, dim);
	
	cout << "Podana macierz:" << endl;
	
	printTridiagonalMatrix(elemSup, elemDia, elemSub, dim);
	
	/*
	 * Sprawdzamy czy macierz jest diag. dominująca.
	 */
	 
	try
	{
		for(int i = 0; i < dim; i++)
		{
			if(i == 0)
			{
				if(!(fabsl(elemDia[i]) > fabsl(elemSup[i])))
				{
					throw SolveException("Macierz nie jest diag. dominujaca!");
				}
			}
			else if (i == dim - 1)
			{
				if(!(fabsl(elemDia[i]) > fabsl(elemSub[i - 1])))
				{
					throw SolveException("Macierz nie jest diag. dominujaca!");
				}
			}
			else
			{
				if(!(fabsl(elemDia[i]) > fabsl(elemSub[i - 1]) + fabsl(elemSup[i])))
				{
					throw SolveException("Macierz nie jest diag. dominujaca!");
				}
			}
		}
		
		/* 
		* Algorytm Thomasa.
		*/
		
		long double m;
		
		int i;
		
		for (i = 1; i < dim; i++)
		{
			m = elemSub[i - 1] / elemDia[i - 1];
			elemDia[i] = elemDia[i] - m * elemSup[i - 1];
			elemB[i] = elemB[i] - m * elemB[i-1];
			
			cout << "Wektor diag. po kroku:" << endl;
			
			printVector(elemDia, dim);
			
			cout << "Wektor wyr. wol. po kroku:" << endl;
			
			printVector(elemB, dim);
			
		}
		
		x[dim - 1] = elemB[dim - 1] / elemDia[dim-1];
		
		for (i = dim - 2; i >= 0; i--)
		{
			x[i] = (elemB[i] - elemSup[i] * x[i + 1]) / elemDia[i];
			
			cout << "Wektor x po kroku:" << endl;
			
			printVector(x, dim);
		
		}
		
		cout << "Wektor x:" << endl;
		
		printVector(x, dim);
		
	}
	catch (SolveException ex)
	{
		cout << ex.getMsg();
		return -1;
	}
	system("PAUSE");
	return 0;
}


// Program implementujący metodę Newtona rozwiązywania układu równań nieliniowych.

#include <iostream>
#include <iomanip>

#include <cmath>

using namespace std;

	/*
	 * Przykładowy układ równań:
	 * 
	 * f1(x, y) = y^3 + 4x = 0
	 * 	 * i f2(x, y) = x^2 + y^2 - 8 = 0
	 * 
	 * Pochodne cząstkowe:
	 * 
	 * d(f1(x, y)) / d(x) = 4
	 * d(f1(x, y)) / d(y) = 3y^2
	 * d(f2(x, y)) / d(x) = 2x
	 * d(f2(x, y)) / d(y) = 2y
	 * 
	 * Wyniki:
	 * 
	 * X1(-2, 2)
	 * X2(2, -2)
	 * 
	 */

long double fun1(long double x, long double y)
{
	return y * y * y + 4.0 * x;
}

long double fun1p1(long double x, long double y)
{
	return 4.0;
}

long double fun1p2(long double x, long double y)
{
	return 3.0 * y * y;
}

long double fun2(long double x, long double y)
{
	return x * x + y * y - 8.;
}

long double fun2p1(long double x, long double y)
{
	return 2.0 * x;
}

long double fun2p2(long double x, long double y)
{
	return 2.0 * y;
}

void printMatrix(long double ** elem, int n, int m)
{
	// cout << setprecision(5) << fixed;
	
	for(int i = 0; i < n; i++)
	{
		cout << "|";
		for(int j = 0; j < m; j++)
		{
			cout << setw(10) << elem[i][j] << " ";
		}
		cout << "|" << endl;
	}
}

void printVector(long double * elem, int n)
{
	// cout << setprecision(5) << fixed;
	
	cout << "[";
	
	for(int i = 0; i < n; i++)
	{
		cout << setw(10) << elem[i] << " ";
	}
	
	cout << "]" << endl;
}

int main(int argc, char* argv[])
{
	long double ** jak;	
	jak = new long double * [2];
	for(int i = 0; i < 2; i++)
	{
		jak[i] = new long double [2];
	}
	
	long double * pocz;
	pocz = new long double[2];
	pocz[0] = 0;
	pocz[1] = sqrtl(8.0);
	
	
	// Równanie:
	// F(X) + Jakobian(X) * H = 0,
	// gdzie F(X) oznacza układ równań,
	// H wektor rozwiązań.
	
	int i = 0; // Numer iteracji.
	int i_max = 100; // Maksymalna liczba iteracji.
	
	long double * stareX;
	stareX = new long double[2];
	stareX[0] = pocz[0];
	stareX[1] = pocz[1];
	
	long double * noweX;
	noweX = new long double [2];
	
	long double * noweH;
	noweH = new long double [2];
	
	long double * noweF;
	noweF = new long double [2];
	
	long double ** elem;
	elem = new long double * [2];
	for(int i = 0; i < 3; i++)
	{
		elem[i] = new long double [3];
	}
	
	long double e = 1e-6;
	
	while((fabsl(fun1(stareX[0], stareX[1]) - fun2(stareX[0], stareX[1])) > e) && i < i_max)
	{
		cout << "Krok " << i + 1 << "." << endl;
		
		// Rozw. układ: Jakobian od stareX, stareY * H = -1 * Wektor,
		// gdzie wsp. są wartości, kolejno, dla funkcji 1 od stareX, stareY,
		// oraz dla funkcji 2 od stareX, stareY,
		
		// Obliczamy wartości jakobianu:
		
		jak[0][0] = fun1p1(stareX[0], stareX[1]);
		jak[0][1] = fun1p2(stareX[0], stareX[1]);
		jak[1][0] = fun2p1(stareX[0], stareX[1]);
		jak[1][1] = fun2p2(stareX[0], stareX[1]);
		
		cout << "Jakobian:" << endl;
		printMatrix(jak, 2, 2);
		
		// Obliczamy wartości funkcji 1 i 2 dla stareX,
		// czyli wyrazy wolne układu równań:
		
		noweF[0] = fun1(stareX[0], stareX[1]);
		noweF[1] = fun2(stareX[0], stareX[1]);
		
		cout << "Wektor wartości funkcji:" << endl;
		printVector(noweF, 2);
		
		//// Tworzymy macierz uzupełnioną o wektor wartości funkcji:
		
		//elem[0][0] = jak[0][0];
		//elem[0][1] = jak[0][1];
		//elem[0][2] = - noweF[0];
		
		//elem[1][0] = jak[1][0];
		//elem[1][1] = jak[1][1];
		//elem[1][2] = - noweF[1];
		
		//cout << "Macierz uzupełniona dla metody Gaussa:" << endl;
		//printMatrix(elem, 2, 3);
		
		//// Metoda Gaussa:
		
		//int dim = 2;
		//long double m;
		//long double e = 1e-6;
		//int i2;
		//int j2;
		//int k2;
		
		//for(i2 = 0; i2 < dim - 1; i2++)
		//{
			//if(fabsl(elem[i2][i2]) < e)
			//{
				//cout << "Błąd!" << endl;
				//break;
			//}
			
			//for(j2 = i2 + 1; j2 < dim; j2++)
			//{
				//m = -elem[j2][i2] / elem[i2][i2];
				
				//for(k2 = 0; k2 < i + 1; k2++)
				//{
					//elem[j2][k2] = 0;
				//}
				//for(k2 = i2 + 1; k2 < dim + 1; k2++)
				//{
					//elem[j2][k2] += m * elem[i2][k2];
				//}
				
				//cout << "Macierz po kroku eleminacji Gaussa:" << endl;
				//printMatrix(elem, dim, dim + 1);
			//}
		//}
	
		//long double * res;
		//res = new long double [dim];
		//long double tmp;
	
		//for(i2 = dim - 1; i2 >= 0; i2--)
		//{
			//if(fabsl(elem[i2][i2]) < e)
			//{
				//cout << "Błąd!" << endl;
				//break;
			//}
			
			//tmp = elem[i2][dim];
			
			//for(j2 = dim - 1; j2 > i2; j2--)
			//{
				//tmp -= elem[i2][j2] * res[j2];
			//}
			//res[i2] = tmp / elem[i2][i2];
		//}
		
		//noweH[0] = res[0];
		//noweH[1] = res[1];
		
		noweH[1] = (jak[0][0] * fun2(stareX[0], stareX[1]) - jak[1][0] * fun1(stareX[0], stareX[1])) / 
		(jak[1][0] * jak[0][1] - jak[0][0] * jak[1][1]);
		
		noweH[0] = -(jak[0][1] * noweH[1] + fun1(stareX[0], stareX[1])) / jak[0][0];
		
		// Podajemy rozwiązanie układu równań liniowych.
		
		cout << "Wektor H:" << endl;
		printVector(noweH, 2);
		
		// Nowe przybliżenia:
		
		for(int j = 0; j < 2; j++)
		{
			noweX[j] = stareX[j] + noweH[j];
		}
		
		cout << "Nowe przybliżenia:" << endl;
		printVector(noweX, 2);
		
		stareX[0] = noweX[0];
		stareX[1] = noweX[1];
		
		i++;
	}
	
	return 0;
}


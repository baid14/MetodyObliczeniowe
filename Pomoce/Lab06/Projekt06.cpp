#include <iostream>
#include <stdio.h>
#include <string.h>
using namespace std;
/******************************************************************************/

//------------------ Wyœwietlanie macierzy -------------------
void wyswietl_macierz(double** macierz, int n, int m)
{
	for(int i=0; i<n; i++)
	{
		cout << "| ";
		for(int j=0; j<m; j++)
		{
			printf(" %.3lf  ",macierz[i][j]);
		}
		cout << " |" << endl;
	}
}

//------------------ Realizacja algorytmu Thomasa -------------------
double *algorytm_thomasa(double **macierz, double *b, int n)
{
	double *beta, *gamma;
	beta = new double[n]; 
    gamma = new double[n];
	double *x = new double[n];

	//------------------ Obliczenie bety i gammy -------------------
	beta[0] = - (macierz[0][1]/macierz[0][0]);
	gamma[0] = (b[0]/macierz[0][0]);

	for(int i=1; i<n; i++)
	{
		if (i<=n-2)
			beta[i] = - ((macierz[i][2+i-1])/((macierz[i][i-1]*beta[i-1])+macierz[i][1+i-1]));

		gamma[i] = (b[i] - (macierz[i][i-1] * gamma[i-1]))/(macierz[i][i-1] * beta[i-1] + macierz[i][1+i-1]);
	}
	beta[n-1] = 0;
    return beta;
}

/******************************************************************************/
/******************************************************************************/
int main()
{
    //------------------ Tworzenie macierzy A oraz nadanei wartoœci -------------------
    double A[5][5] = {10,4,0,0,0,
                      1,11,3,0,0,
                      0,2,12,2,0,
                      0,0,3,13,1,
                      0,0,0,4,14};
                      
    //------------------ Tworzenie wektora b oraz nadanei wartoœci -------------------
    double b[5] = {18,32,44,36,22};
    
    double** tabA = new double *[5];
    
    for(int i=0; i<5; i++)
    {
      tabA[i] = new double[5];
      memcpy(tabA[i], A[i], sizeof(double)*5);
    }
    
    double* tabB = new double [5];           
	memcpy(tabB,b,sizeof(double)*5);

	double* x = new double [5];
	
	//------------------ Wyœwietlenie macierzy A -------------------
	cout<<"Macierz trojprzekatniowa \n"<<endl;
	wyswietl_macierz(tabA, 5, 5);
	
	//------------------ Wyœwietlenie wektora b -------------------
	cout<<"\nWektor b \n"<<endl;
	
 	for(int i = 0; i<5; i++)
		cout << "| " << tabB[i] <<  " |" <<endl;
		
    //------------------ Zastosowanie algorytmu Thomasa -------------------
    x = algorytm_thomasa(tabA, tabB, 5);
    
    //------------------ Wyœwietlenie wektora x - Rozwi¹zanie uk³adu Ax = b za pomoc¹ alg. Thomasa -------------------
    cout<<"\nRozwiazanie uk³adu Ax=b za pomoc¹ metody Thomasa"<<endl;
	cout<<"--- Wektor x \n"<<endl;
	for(int i = 0; i<5; i++)
		cout << "| " << x[i] << "  |"<< endl;
                      
}

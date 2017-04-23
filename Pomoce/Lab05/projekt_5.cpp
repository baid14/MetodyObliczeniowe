#include <iostream>
#include <iomanip>
#include <math.h>
using namespace std;


void wypiszMac(double **mac, int n)
{
	for(int i=0; i<n; i++)
	{
		cout << "| ";
		for(int j=0; j<n; j++)
		{
			cout << setw(5) << mac[i][j];
		}
		cout << " |" << endl;
	}
}

void wypiszWek(double *wek, int n)
{
	for(int i=0; i<n; i++)
	{
		cout << "| " << setw(3) << wek[i] << " |" << endl;
	}
}

double **macierz(int n)
{
	double **macierz = new double*[n];
	for(int i=0; i<n; i++)
		macierz[i] = new double[n];
	return macierz;
}


long double wyznacznik(int st, int wiersz, int *wek, double **mac)
{  
  int kolumny[10]; // wektor kolumn dla podmacierzy
  long double wynik;
  int i,j,k,m;

  if(st == 1)
    return mac[wiersz][wek[0]];
  else
  {
    wynik = 0; m = 1;
    for(i = 0; i < st; i++)
    {
      k = 0;
      for(j = 0; j < st - 1; j++)
      {
        if(k == i) k++;
        kolumny[j] = wek[k++];
      }
      wynik += m * mac[wiersz][wek[i]] * wyznacznik(st - 1, wiersz + 1, kolumny, mac);
      m = -m;
    }
  return wynik;
  }
}

void pods_ele(double **A, double *b, int i, int n)
{
	double eps=0.000001;
	double* p_a;  
	double p_b;
	double max = fabs(A[i][i]);
	int max_tab = i;
	int j; 
	for(j=i; j<n; ++j)
	      {
               if( fabs(A[j][j]) > max)
                   {
                   max = fabs(A[j][j]);
                   max_tab = j;
                   }
          }
    if( max < eps ) // sprawdzam czy to macierz osobliwa, jeœli tak to det = 0
        {
        printf("Macierz osobliwa!");
        system("PAUSE");
        exit(0);    
        }
    p_a = A[max_tab]; 
    A[max_tab] = A[i];
    A[i] = p_a;
    
	p_b = b[max_tab];
    b[max_tab] = b[i];
    b[i] = p_b;
}

void rozklad_LU(double **A, double* b, int n)
{
     double eps = 0.0001;
	double x;
	for(int k=0; k<n-1; k++)
	{
            
		for(int i=k+1; i<n; i++)
		{
            if( A[k][k] < eps )
                pods_ele( A,b, k, n );
			x = A[i][k]/A[k][k]; // dzielimy elementy przez element na przek¹tnej
			A[i][k] = x;
			for(int j=k+1; j<n; j++)
			{
				A[i][j] = A[i][j] - (x*A[k][j]);
			}
		}
	}
}

//z to Wektor rowi¹zañ uk³adu równan Lz = b

void rozwiaz_uklad(double **A, double *b, int n)
{
	double *z, *x;
	z = new double[n];
	x = new double[n];
	double suma;

    //obliczanie wektora z
	for(int i=0; i<n; i++)
	{
		suma = 0;
		for(int j=0; j<=i-1; j++)
		{
			suma += A[i][j]*z[j];
		}

		z[i] = b[i]-suma;
	}


	 cout << "\n\nLz=b.\nWektor z:"<<endl;
     wypiszWek(z, n);

	 //obliczanie rozwi¹zañ
	for(int i=n-1; i>=0; i--)
	{
		suma = 0;
		for(int j=i+1; j<n; j++)
		{
			suma +=A[i][j]*x[j];
		}

		x[i] = (z[i]-suma)/A[i][i];
	}

	 cout << "\n\nUx=z.\nWektor x:"<<endl;
     wypiszWek(x, n);
}

int main()
{
	double **A, *b;
	int n = 4;				
	A = macierz(n);
	b = new double[n];
     
     A[0][0]=1.0; 
	 A[0][1]=2.0; 
	 A[0][2]=2.0; 
	 A[0][3]=1.0;
     A[1][0]=2.0; 
	 A[1][1]=4.0; 
	 A[1][2]=4.0; 
	 A[1][3]=1.0;
     A[2][0]=2.0; 
	 A[2][1]=2.0; 
	 A[2][2]=2.0; 
	 A[2][3]=1.0;
     A[3][0]=1.0; 
	 A[3][1]=1.0; 
	 A[3][2]=2.0; 
	 A[3][3]=1.0;
    
	 b[0]=1.0;
	 b[1]=2.0;
	 b[2]=3.0;
	 b[3]=4.0; 
     
	cout << "Macierz A" << endl;
	wypiszMac(A, n);

	cout << "\nWektor b" << endl;
	wypiszWek(b,n);

	int *tab = new int[n];
    for(int i=0; i<n; i++)
    {
            tab[i] = i;
    }

    if( fabs( wyznacznik(n,0,tab,A )) < 0.0000001) 
    {
            cout << "Podana macierz jest osobliwa" << endl;
            exit(1);
    }
	rozklad_LU(A, b, n);
	cout << "Macierz LU" << endl;
	wypiszMac(A, n);
	rozwiaz_uklad(A,b,n);
	system("PAUSE");
}




#include <iostream>
#include <iomanip>
#include <math.h>
using namespace std;
/******************************************************************************/

//------------------ Tworzenie wektora -------------------
double *nowy_wektor(int n)
{
       double *wektor = new double [n];
       return wektor;
}

//------------------ Tworzenie macierzy--------------------
double **nowa_macierz(int n)
{
       double **macierz = new double *[n];
       
       for(int i=0; i<n; i++)
          macierz[i] = new double[n];
       
       return macierz;
}

//------------------ Wyœwietlanie wektora -------------------
void wyswietl_wektor(double *b, int n)
{
     for(int i=0; i<n; i++)
     {
       cout<<setw(4)<<b[i]<<endl;
     }
}

//------------------ Wyœwietlanie macierzy -------------------
void wyswietl_macierz(double **A, int n)
{
     for(int i=0; i<n; i++)
     {
        for(int j=0; j<n; j++)
        {
          cout<<setw(5)<<A[i][j];
        }
        cout<<endl;
     }
}

//------------------ Wartoœæ wyznacznika -------------------
long double wyznacznik(int stopien, int wiersz, int *wk, double **A)
{
     int m, k;
     int kolumny[10];                   //Wektor kolumn dla podmacierzy
     long double wynik;
     
     if(stopien == 1)
       return A[wiersz][wk[0]];
     else
     {
         wynik = 0;
         m = 1;
         
         for(int i=0; i<stopien; i++)
         {
           k=0;
           
           for(int j=0; j<stopien-1; j++)
           {
             if(k == i)
               k++;
             kolumny[j] = wk[k++];   
           }
           wynik = wynik+m*A[wiersz][wk[i]]*wyznacznik(stopien-1, wiersz+1, kolumny, A);
           m=-m;
         }
         
         return wynik;
     }
}

//------------------ Wybór czêœciowego elementu podstawowego -------------------
void element_podstawowy(double **A, double *b, int i, int n)
{
     double eps = 0.000001;
     
     double* swap;
     double swap_b;
     double max = fabs(A[i][i]);
     int max_tab = i;
     
     for(int j=i; j<n; ++j)
     {
       if(fabs(A[j][j]) > max)
       {
         max = fabs(A[j][j]);
         max_tab = j;
       }
     }
     
     if(max < eps)
     {
       printf("Macierz osobliwa!");
       system("pause");
       exit(0);
     }
     
     swap = A[max_tab];
     swap_b = b[max_tab];
     A[max_tab] = A[i];
     b[max_tab] = b[i];
     A[i] = swap;
     b[i] = swap_b;
     
}

//------------------ Rozk³ad LU -------------------
void rozklad_LU(double **A, double *b, int n)
{
  double eps = 0.0001;
  double x;
  
  for(int k=0; k<n-1; k++)
  {
   for(int i=k+1; i<n; i++)
   {
    if(A[k][k] < eps)
      element_podstawowy(A, b, k, n);
    
    x = A[i][k]/A[k][k];
    A[i][k]=x;
    
    for(int j=k+1; j<n; j++)
    {
      A[i][j] = A[i][j] - (x*A[k][j]);
    }
   }
  }
}

//------------------ Rozwiazanie uk³adu równañ Ly=b, gdzie y jest wektorem rozwiazañ -------------------
void rozwiaz(double **A, double *b, int n)
{
     double *y, *x;
     double suma;
     
     y = nowy_wektor(n);
     x = nowy_wektor(n);
     
     // --- WEKTOR y ----
     for(int i=0; i<n; i++)
     {
       suma = 0;
       for(int j=0; j<=i-1; j++)
       {
         suma = suma + A[i][j]*y[j];
       }
       
       y[i] = b[i]-suma;
     }
     
     cout << "Rozwiazaniem ukladu Ly=b jest\nWektor y"<<endl;
     wyswietl_wektor(y, n);
     
     // --- WEKTOR x ----
     for(int i=n-1; i>=0; i--)
     {
       suma = 0;
       for(int j=i+1; j<n; j++)
       {
         suma = suma + A[i][j]*x[j];
       }
       x[i] = (y[i] -suma)/A[i][i];
     }
     
     cout << "Rozwiazaniem ukladu Ux=y jest\nWektor x"<<endl;
     wyswietl_wektor(x, n);
}

/******************************************************************************/
/******************************************************************************/
int main()
{
    double *b;                           // Wektor.                                                          
    double **A;                          // Macierz.
    
    int n = 4;                           //Rozmiar.
    
    //------------------ Tworzenie wektora i macierzy -------------------
    A=nowa_macierz(n);
    b=nowy_wektor(n);
    
    //------------------ Nadanie wartoœci dla wektora b -------------------
    b[0]=1.0;   
    b[1]=2.0;   
    b[2]=3.0;   
    b[3]=4.0;
    
    //------------------ Nadanie wartoœci dla macierzy A -------------------
    A[0][0]=1.0;  A[0][1]=2.0;  A[0][2]=2.0;  A[0][3]=1.0;
    A[1][0]=2.0;  A[1][1]=4.0;  A[1][2]=4.0;  A[1][3]=1.0;
    A[2][0]=2.0;  A[2][1]=2.0;  A[2][2]=2.0;  A[2][3]=1.0;
    A[3][0]=1.0;  A[3][1]=1.0;  A[3][2]=2.0;  A[3][3]=1.0;
    
    //------------------ Wyœwietlenie wektora b -------------------
    cout << "Wektor b"<<endl;
    wyswietl_wektor(b, n);
    
    //------------------ Wyœwietlenie macierzy A -------------------
    cout<<"\nMacierz A"<<endl;
    wyswietl_macierz(A, n);
    
    
    //------------------ Tworze nowy wektor -------------------
    int *tab = new int[n];
    
    for(int i=0; i<n; i++)
    {
      tab[i] = i;
    }
    
    if( fabs( wyznacznik(n, 0, tab, A)) < 0.0000001)
    {
        cout<< "Podana macierz jest osobliwa"<<endl;
        exit(1);
    }
    
    //------------------ Stosuje rozk³ad LU-------------------
    rozklad_LU(A, b, n);
    cout<<"\nMacierz LU"<<endl;
    wyswietl_macierz(A, n);
    
    //------------------ Rozwi¹zanie uk³adu-------------------
    rozwiaz(A, b, n);
    
    system("pause");
}

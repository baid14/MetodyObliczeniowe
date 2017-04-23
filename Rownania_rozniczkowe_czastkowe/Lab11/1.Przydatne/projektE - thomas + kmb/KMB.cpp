#include<iostream>
#include<cstdio>
#include<cmath>
#include<cstdlib>
#include <iomanip>
#include<fstream>
#include "calerf.h"
using namespace std;

//====================================== Sta³e ========================================
const double D=1.0;
const double tmax=2.0;
const double a=6*sqrt(D*tmax);
const double lambda=0.4;
const double h=0.1;     //kluczowa sta³a
const double dt=(lambda*(h*h)/D);
const int N = static_cast<const int>((tmax/dt)+1);
const int M = static_cast<const int>((a*2/h)+1); 

//============================ Program korzysta z funkcji =============================
double *nowy_wektor(int n);
double **nowa_macierz( int n, int m ); 
void delete_matrix(double **matrix, int n);
void delete_vector(double *x);
void print_vector(double *x, int k);
void print_matrix(double **x, int r, int c); 
void fill_dt(double *x, double delta, int k);
void fill_h(double *x, double delta, int k);
void save_as_file( char* filename, double **matrix, int r, int c);
double analityczne(double **ANAL, int N, int M, double *dt_levels, double *h_levels, double dt, double h);
void klas_met_bezposred(double **KMB, int N, int M, double *dt_levels, double *h_levels, double dt, double h, double lambda);
void blad(double **BLAD, double **ANAL, double **KMB, int N, int M);

//=====================================================================================
int main()
{
    cout.setf( ios::fixed, ios::floatfield );
    cout.precision(1);
    
    
    double *dt_levels = nowy_wektor(N);
    double *h_levels = nowy_wektor(M);
    fill_dt(dt_levels, dt, N);
    fill_h(h_levels, h, M);
     
    // zapis poziomow czasowych i przestrzennych do pliku:
    char *filename = "hdt_KMB.txt";
    ofstream file( filename );
    if( !file.is_open() )
    {
        cout << "Nie moge zapisac: " << filename <<" !!!"<< endl;
        exit(1);
    }
    file.setf(ios::fixed, ios::floatfield); 
    file.precision(4); 
    file << endl << "      POZIOMY dt i h:             " << endl << endl;
    file << "dt\\h: \t";
    for( int i=0; i<M; i++)
         file  << h_levels[i]<< "\t";
    file << endl;
    for(int i=0; i<N; i++)
         file <<  dt_levels[i]<<endl;
    file.close();
    cout << "Zapisano: " << filename  << endl;
    
    //--------------------------------------------------------------------------------------
       
    double **ANAL = nowa_macierz(N, M);
    analityczne(ANAL, N, M, dt_levels, h_levels, dt, h);
    
    double **KMB = nowa_macierz( N, M );
    klas_met_bezposred(KMB, N, M, dt_levels, h_levels, dt, h, lambda);
    
    double **BLAD = nowa_macierz(N,M);
    blad(BLAD, ANAL, KMB, N, M);
    
    delete_matrix(KMB, N);
    delete_matrix(ANAL, N);
    delete_vector(dt_levels);
    delete_vector(h_levels);
    
    system("PAUSE");
return 0;   
}
double *nowy_wektor(int n)
{ 
       double *x = new double[n]; 
       return x; 
}

double **nowa_macierz( int n, int m ) // n - liczba wierszy, m - liczba kolumn
{
       double **matrix=new double *[n];
       if(!matrix) 
       { 
           cout << "!Blad alokacji!"; 
           exit(1);
       } 
       for(int i=0; i<n; i++) 
       { 
           matrix[i] = new double [m]; 
           if( !matrix[i] ) 
           { 
               cout << "!Blad alokacji!"; 
               exit(1); 
           } 
       } 
       return matrix;
}

void delete_matrix(double **matrix, int n)
{
     for(int i=n-1; i>0; i--)
          delete []matrix[i];
     delete []matrix;
}

void delete_vector(double *x)
{
     delete []x;
}

void print_vector(double *x, int k)
{
     for( int i = 0; i < k; i++ )
          cout << x[i] << " | ";
     cout << endl;
}

void print_matrix(double **x, int r, int c) // matrix[c][r]
{
     for(int i=0; i<r; i++)
     {
          for(int j=0; j<c; j++)
               cout << x[i][j] << "    ";
          cout << endl;
     }
}

void fill_dt(double *x, double delta, int k)
{
     x[0] = 0; // wynika z zalozen
     for(int i=1; i<k; i++)
          x[i]=x[i-1]+delta;
}

void fill_h(double *x, double delta, int k)
{
     x[0] = -a; // wynika z warunkow poczatkowych
     for(int i=1; i<k; i++)
          x[i]=x[i-1]+delta;
}

void save_as_file( char* filename, double **matrix, int r, int c)
{
     ofstream file( filename );
     if( !file.is_open() )
     {
         cout << "Nie moge zapisac: " << filename <<" !!!"<< endl;
         exit(1);
     }
     file.setf(ios::scientific, ios::floatfield); 
     file.precision(4); 
     file << endl << "                " << filename << ":                " << endl << endl;
     for( int i = 0; i < r; i++ )
     {
          file << endl;
          for( int j = 0; j < c; j++ )
              file << matrix[i][j]<<"\t";
     }
     file.close();
     cout << "Zapisano: " << filename  << endl;
}

// Rozwiazanie analityczne U(x,t): 
double analityczne(double **ANAL, int N, int M, double *dt_levels, double *h_levels, double dt, double h)
{ 
     for(int i=0; i<N; i++)    //warunek brzegowy U(-INF,t)=1
     ANAL[i][0] = 1.;
     for(int i=0; i<N; i++)    //warunek brzegowy U(+INF,t)=0
     ANAL[i][M-1]=0.;
     for(int i=0; i<M; i++)  //warunek poczatkowy  U(x,0)= (1 dla x<0) (1 dla x>=0) 
        if(h_levels[i]<0.) ANAL[0][i]=1.;
                    else ANAL[0][i]=0.;
      
       for(int i=1; i<N; i++)
       {
            for( int j=1; j<M; j++) 
                 ANAL[i][j]=0.5*erfc((h_levels[j])/(2*sqrt(D*dt_levels[i])));
       }
       save_as_file("an_KMB.txt", ANAL, N, M);
}

void klas_met_bezposred(double **KMB, int N, int M, double *dt_levels, double *h_levels, double dt, double h, double lambda)
{
     for(int i=0; i<N; i++)    //warunek brzegowy U(0,t)=0
     KMB[i][0] = 1;
     for(int i=0; i<N; i++)    //warunek brzegowy U(INF,t)=0
     KMB[i][M-1]=0;
     for(int i=0; i<M; i++)  //warunek poczatkowy  U(x,0)= (1 dla x<0) (1 dla x>=0) 
        if(h_levels[i]<0.) KMB[0][i]=1.;
                    else KMB[0][i]=0.;
     
     for(int k = 1; k < N; k++){
        for(int i = 1; i < M-1; i++){
             KMB[k][i]=KMB[k-1][i]+lambda*(KMB[k-1][i-1]-(2*KMB[k-1][i])+KMB[k-1][i+1]);
        }
     }
     save_as_file("num_KMB.txt", KMB, N, M);
}

void blad(double **BLAD, double **ANAL, double **KMB, int N, int M)
{
    for(int i=0; i<N; i++)
    {
            for(int j=0; j<M; j++) 
            {
                BLAD[i][j]=fabs(KMB[i][j]-ANAL[i][j]);
            }    
            }
            save_as_file("bledy_KMB.txt", BLAD, N, M);
}

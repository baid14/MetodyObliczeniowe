#include <iostream>
#include <stdlib.h>
#include <cstdlib>
#include <fstream>
#include <String>
#include <Math.h>
#include <sys/timeb.h>
using namespace std;

fstream file,file2,file3,file4,file5,file6,file7,file8,file9,file10;
	


//deklaracja sta³ych:
const double D = 1.0;
const double lambda_bez = 0.4; //dla metod bezposredniej
const double lam_pos = 1.0; //dla metod posrednich
const double h = 0.1; //krok dla x
const double dt = (lambda_bez*h*h)/D; //krok dla t
const double PI=3.1415;
const double x_max=1; //maksymalny przedzial x
const double x_min=0; //minimalny przedzial x
const double t_max=0.5; //maksymalny przedzial x
const double t_min=0; //minimalny przedzial x


//wymiary macierzy do obliczeñ numerycznych:
//const int N = (int)((2*x_max/h)+1);//liczba kolumn (os x)
const int N = 40; //liczba kolumn(x)
const int M = (int)((t_max/dt)+1);//liczba wierszy (os t)

//const int M=(t_max-t_min)/dt;
//const int N=(x_max-x_min)/h;

//**************************************************************************************************
//tworzy now¹ macierz MxN
double **macierz(int m, int n) 
{
       double **a;
       a=new double *[m];
       for(int i=0; i<m; i++)
       a[i]=new double[n];
       return a;
}

//**************************************************************************************************
//usuwa macierz
void usun_macierz(double **a, int m)
{
     for(int i=m-1; i>=0; i--)
     delete []a[i];
     delete []a;
}

//**************************************************************************************************
//zapisywanie macierzy do pliku
void zapisz(double **macierz, fstream &plik)
{
     
     //plik.setf(ios::scientific, ios::floatfield); 
     //plik.precision(4); 
     plik << ";";
     for(int i = 0; i<N; i++){
	  		 plik << i*h << ";";
	  	}
	  	plik << endl;
             for(int j = 0; j<M; j++){
             	plik << j*dt << ";";
				for(int i=0;i<N;i++){
					plik << macierz[j][i] << ";";
				}
				plik<<"\n";
			 }
                
     plik.close();
}

//**************************************************************************************************
void zapiszW(double *wek, fstream &plik)
{ 
             for(int j = 0; j<M; j++){
					plik << j*h << ";" << wek[j] << "\n";
				}             
     plik.close();
}

/*********************************************************************************************/
/*metoda obliczaj¹ca czas wykonywania obliczen*/
double Czas()
{
       struct timeb czas;
       double sekundy;
       ftime(&czas);

       sekundy= (double) czas.time;
       sekundy += (double) czas.millitm / 1000.0;
       return sekundy;
}

//**************************************************************************************************
//Metoda rozwi¹zuj¹ca analitycznie równanie rózniczkowe
void analitycznie(double **A, double *x, double *t)
{
	//warunek pocz¹tkowy
     for(int i=0;i<N;i++){
             A[0][i]=sin(PI*x[i]); 
     }
	 
	 //warunek brzegowy
     for(int i=0;i<M;i++){
             A[i][N-1]=0; 
             A[i][0]=0; 
	}
	
	//wype³nienie ca³ej macierzy:
     for(int i =1; i<M; i++){
             for(int j =1; j<N-1; j++){
                     A[i][j] = exp(-(PI*PI)*D*t[i])*sin(PI*x[j]);
					 }
        }
     
     //zapisanie rozwi¹zania do pliku:
     zapisz(A, file);
}

//**************************************************************************************************
//klasyczna metoda bezposrednia
void klas_bezp(double **KMB, double *x, double *t, double lambda_bez)
{
 	 file8.open("Czas_KMB.txt",ios::out);
 	 double start, stop;
     start = Czas();
	  //warunek pocz¹tkowy
     for(int i=0;i<N;i++){
             KMB[0][i]=sin(PI*x[i]); 
     }
	 
	 //warunek brzegowy:
     for(int i =0; i<M; i++){
             KMB[i][N-1] = 0;
			 KMB[i][0] = 0;
             
     }   
     //wype³nienie ca³ej macierzy:
     for(int k = 1; k < M; k++){
        for(int i = 1; i < N-1; i++){
             KMB[k][i]=KMB[k-1][i]+lambda_bez*(KMB[k-1][i-1]-(2*KMB[k-1][i])+KMB[k-1][i+1]);
        }
     }
	 stop = Czas();
	 file8 << "Czas[sek]: \t" << (stop - start )<< "\n";
	 file8 << "H: \t"<< h << "\n";
	 file8 << "N: \t"<< N << "\n";	 
     //zapisanie rozwi¹zania do pliku:
     zapisz(KMB, file2);
}

//**************************************************************************************************
void blad(double *BLAD, double **KMB, double **A, fstream &plik)
{
 	 double temp;
     for(int i =0; i<M; i++)
     {
	  		 double max = 0.0;
             for(int j =0; j<N; j++){
                     temp = fabs(KMB[i][j] - A[i][j]);
                     if(temp>max){
					 			  max = temp;
					 			  }
                     }
                     BLAD[i] = max;
     }

     zapiszW(BLAD, plik);
     
}

//**************************************************************************************************
double norma(double *dane){
	double w = dane[0];
	for( register int i = 1; i < N; i++ ){
		if( w < dane[i]  )
			w = dane[i] ;
	}
	return w;
}

//**************************************************************************************************
//funkcja dokonuj¹ca czêœciowego wyboru elementów podstawowych
void szukaj_zamien(double **A, int *indeksy, int start){
	int max = start;     
	for(register  int i = start+1; i < M; i++ )
		if( fabs( A[indeksy[i]][start] ) > A[indeksy[max]][start] )
			max = i;  
	//sprawdzanie czy uda³o siê znaleœæ odpowiedni element          
	if( A[indeksy[max]][start] == 0 ){
		printf("Bladna macierz A\n");
		exit(1);    
	}             
	indeksy[start] = max;
	indeksy[max] = start;
}

//**************************************************************************************************
//funkcja dokonuj¹ca dekompozycji LU
void rozklad_LU(double **LU, double **A, int *indeksy){
	double wsp;
	for(int i = 0; i < N; i++ ){
		for(int j = 0; j < N; j++ ){
			LU[indeksy[i]][j] = A[indeksy[i]][j];
		}
	}
	for( int i = 0; i < N; i++ ){  
		for(int j = i+1; j < N; j++ ){
			if(LU[indeksy[i]][i] == 0 ) szukaj_zamien(LU,indeksy,i);
			wsp = LU[indeksy[j]][i]/LU[indeksy[i]][i];
			for(int k = i; k < M; k++ ){
				LU[indeksy[j]][k] = LU[indeksy[j]][k] - LU[indeksy[i]][k]*wsp;    
				if( i < j ){
					LU[indeksy[j]][i] = wsp;      
				}
			}     
		}
	}
}

//**************************************************************************************************
//funkcja dokonuj¹ca obliczeñ wektora y z równania Ly=b
void obliczanie_y(double **LU, double *Y, double *B, int *indeksy, const int n){
	double S;
	for(register  int i = 0; i < n; i++ ){
		S = 0;
		for(register  int j = 0; j < i; j++ ){
			S += Y[j]*LU[indeksy[i]][j];
		}
		Y[i] = B[indeksy[i]] - S;   
	}     
}

//**************************************************************************************************
//funkcja dokonuj¹ca obliczeñ wektora x z równania Ux=y
void obliczanie_x(double **LU, double *X, double *Y, int *indeksy, const int n){
	double S;
	register int j;
	for( register int i = n-1; i >= 0; i-- ){
		S = 0;
		for( j = n-1; j > i; j-- )
			S += X[j]*LU[indeksy[i]][j];
		X[i] = (Y[i] - S)/LU[indeksy[i]][i];   
	}     
}

//**************************************************************************************************
//funkcja inicjalizuj¹ca tablicê indeksów potrzebna dla powy¿szych funkcji.
void init_kolumny(int *indeksy, const int n ){
	for( register int i = 0; i < n; i++ )
		indeksy[i] = i;     
}

//**************************************************************************************************
void crank_nicolson_LU(double **A,double *x, double *t, double lam_pos){
	register int i,j,k;
	double *b,*z;

	double **LU = new double * [N];
	double **WSP = new double * [N];
	for( i = 0 ; i < N ; i++ ){
		LU[i] = new double [N];
		WSP[i] = new double [N];
	}
	double *Y = new double [N];
	int *indeksy = new int [N];

	b=new double[N];
	z=new double[N];

	for(j=0;j<N;j++)
		A[0][j]=sin(PI*x[i]);

	for(i=0;i<M;i++){
		A[i][0]=0;
		A[i][N-1]=0;
	}

	for( i = 0; i < N; i++ )
		for(  j =0; j< N; j++ ){
			if( j == i-1 || j == i+1){
				WSP[i][j] = -lam_pos/2.0;
			}else if( j == i ){
				WSP[i][j] = 1 + lam_pos;
			}else{
				WSP[i][j] = 0;
			}
		}
	
		init_kolumny(indeksy,N);
		rozklad_LU(LU,WSP,indeksy);
		for( i=0;i<M-1;i++){
			b[0]=b[N-1]=0;
			for(j=1;j<N-1;j++)
				//b[j]=(lambda/2.0)*A[i][j-1]+(1-lambda)*A[i][j]+(lambda/2.0)*A[i][j+1];
				b[j]=-((lam_pos/2.0)*A[i][j-1]+(1.0-lam_pos)*A[i][j]+(lam_pos/2.0)*A[i][j+1]);
			//for(k=0;k<m;k++)
				//b[k]=b[k]+D*dt;

	
			obliczanie_y( LU,Y, b,indeksy,N);
			obliczanie_x( LU,z, Y,indeksy,N);

			for(j=1;j<N-1;j++)
				A[i+1][j]=z[j];
		}
		
		zapisz(A, file4);
}
//**************************************************************************************************
//gauss Seidel
/*void GS(double **A, double *x, double *b){
	double tmp;
	double *xp = new double[N];
	double *tmp_w = new double[N];
	for( int i = 0; i < N; i++ )
			xp[i] =x[i] + 1;
	do{
		
		for( int i = 0; i < N; i++ ){
			tmp = 0;
			for( int j = 0; j < i-1; j++ ){
				tmp += A[j][i] * x[j];
			}
			for( int j = i+1; j < N; j++ ){
				tmp += A[j][i] * x[j];
			}
			x[i] = b[i] - tmp;
			x[i] /= A[i][i];
		}
		for( int i = 0; i < N; i++ ){
			tmp_w[i] = fabs(x[i] - xp[i]);
			xp[i] = x[i];
		}
	}while( norma(tmp_w) > 0.00000000001  );
}


void crank_nicolson_GS(double **A, double *x, double *t, double lam_pos){
	register int i,j,k;
	double *b,*z;

	double **LU = new double * [N];
	double **WSP = new double * [N];
	for( i = 0 ; i < N ; i++ ){
		LU[i] = new double [N];
		WSP[i] = new double [N];
	}

	b=new double[N];
	z=new double[N];


	for(j=0;j<N;j++){
		A[0][j]=sin(PI*x[i]);
		z[j] = 1;
	}

	for(i=0;i<M;i++){
		A[i][0]=0;
		A[i][N-1]=0;
	}

	for( i = 0; i < N; i++ )
		for(  j =0; j< N; j++ ){
			if( j == i-1 || j == i+1){
				WSP[i][j] = -lam_pos/2.0;
			}else if( j == i ){
				WSP[i][j] = 1 + lam_pos;
			}else{
				WSP[i][j] = 0;
			}
		}
		for( i=0;i<M-1;i++){
			b[0]=b[N-1]=0;
			for(j=1;j<N-1;j++)
				//b[j]=(lambda/2.0)*A[i][j-1]+(1-lambda)*A[i][j]+(lambda/2.0)*A[i][j+1];
				b[j]=(lam_pos/2.0)*A[i][j-1]-(lam_pos-1)*A[i][j]+(lam_pos/2.0)*A[i][j+1];
			//for(k=0;k<m;k++)
				//b[k]=b[k]+D*dt;

			GS(WSP,z,b);

			for(j=1;j<N-1;j++)
				A[i+1][j]=z[j];
		}
		zapisz(A, file6 );
}
*/


//**************************************************************************************************
int main(){
	file.open("Analit.csv",ios::out);
		file2.open("BMK.csv",ios::out);
		file3.open("Bledy_BMK.csv",ios::out);
	file6.open("CN-LU.csv",ios::out);
		file5.open("Bledy_CN-LU.csv",ios::out);
//	file4.open("CN-GS.csv",ios::out);
//	file7.open("Bledy_CN-GS.csv",ios::out);
	

	double *x = new double[N];
    double *t = new double[M];
    
    
    for(int i =1; i<N; i++)
            x[i] = x_min + (i*h);
    
    for(int i =1; i<M; i++)
            t[i] = t_min + (i*dt);

	double **A = macierz(M, N); //macierz do rozwi¹zania analitycznego
	double **KMB = macierz(M,N);
	double *BLAD = new double[M];
	double **CN_LU = macierz(M,N);
//	double **CN_GS = macierz(M,N);
	analitycznie(A,x,t);
	klas_bezp(KMB,x,t,lambda_bez);
	blad(BLAD,KMB,A,file3);
//	crank_nicolson_GS(CN_GS,x,t,lambda_bez);
//	blad(BLAD,CN_GS,A,file7);
	//crank_nicolson_LU(CN_LU,x,t,lambda_bez);
	//blad(BLAD,CN_LU,A,file5);
	
	return 0;
}

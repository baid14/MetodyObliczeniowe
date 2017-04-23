#include <iostream>
#include <iomanip>
#include <string.h>
#include <stdio.h>
#include <malloc.h>
#include <math.h>
using namespace std;
/******************************************************************************/
class ukladRownan
{
    //------------------ Sk³adowe -------------------
   static double przyblizenie_poczatkowe[2];
    static int liczba_iteracji;		
    static const int n = 2;	
	double **A;	
	double *b;
    enum kryt { XN_MIN_XPOP = 1 ,MAX_IT = 2, RESIDUUM = 4};
    static short int kryterium;
    double static tolf;                   //Tolerancja f
    double static tolx;                   //Tolerancja x
    
	public:
         
    //------------------ Konstruktor domyœlny -------------------       
	ukladRownan()
	{
        A = new double*[2];
        for( int i = 0; i < 2; ++i )
             A[i] = new double[2];
             
        b = new double[2];
		A[0][0] = 5; A[0][1] = 1;
		A[1][0] = 4; A[1][1] = 10;

		b[0] = 49; 
		b[1] = 30;
    }
    
    //------------------ Metody -------------------
    void drukuj();
    void metodaJacobiego();
    void metodaGaussaSeidela();
    void metodaSOR(double parametr = 0.5);
    double residuum( double x[2] );
    double warunek();
};
//------------------ Dane skladowych -------------------
    double ukladRownan::przyblizenie_poczatkowe[2] = {1,1};        
    int ukladRownan::liczba_iteracji = 100;
    double ukladRownan::tolf = 0.00000001;
    double ukladRownan::tolx = 0.00000001;
    short int ukladRownan::kryterium = RESIDUUM;   //XN_MIN_XPOP = 1 ,MAX_IT = 2, RESIDUUM = 4
    
//------------------ Metoda drukuj¹ca -------------------
void ukladRownan::drukuj()
{
	cout <<"Macierz A:" << endl;
	for(int i=0; i<2; i++)
	{
		cout << "|";
		for(int j=0;j<2; j++)
		{
			cout << setw(3) << A[i][j];
		}
		cout << "  |" << endl;
	}

	cout <<"\nWektor b:" << endl;
	for(int i=0; i<2; i++)
	{
		cout << "| " << setw(3) << b[i] << " |" << endl;
	}
}

//------------------ Metoda JACOBIEGO -------------------
void ukladRownan::metodaJacobiego()
{
 double wynik[2];
 double poprzedni[2] = {przyblizenie_poczatkowe[0], przyblizenie_poczatkowe[1]};
 int iter = 0;
 
 cout<<"Metoda JACOBIEGO\n"<<endl;
	while(1)
	{
     ++iter;
		for(int j=0; j<2; j++)
		{
          wynik[j] = 0;

          for(int k=0; k<2; k++)
          {
            if(k!=j)
             wynik[j] += (-1)*(A[j][k]/A[j][j])*poprzedni[k];
          }
          
          wynik[j] += b[j]/A[j][j];
         }
     
     if( kryterium & XN_MIN_XPOP && max( fabs( wynik[0] - poprzedni[0]), fabs( wynik[1] - poprzedni[1]) ) < tolx )
       break;
     if ( kryterium & MAX_IT && iter >= liczba_iteracji )
       break;
     if( kryterium & RESIDUUM && fabs(residuum(wynik)) < tolf )
       break;
       
     memcpy(poprzedni,wynik,sizeof(double)*2);
	 printf("Iteracja %d\t x1 = %g\t x2 = %g\t", iter, wynik[0], wynik[1] );
     printf("Residuum: %g\n", residuum( wynik ) );
	}
}

//------------------ Metoda GAUSSA-SEIDELA -------------------
void ukladRownan::metodaGaussaSeidela()
{
 double wynik[2];
 double poprzedni[2] = {przyblizenie_poczatkowe[0], przyblizenie_poczatkowe[1]};
 int iter = 0;
 
 cout<<"Metoda GAUSSA-SEIDELA\n"<<endl;
	//wyliczanie przyblizen rozwiazan
	while(1)
	{
     ++iter;
		for(int j=0; j<n; j++)
		{
          wynik[j] = 0;
          
          for(int k=0; k<n; k++)
		  {
            if(k!=j)
             wynik[j] += (-1)*(A[j][k]/A[j][j])*poprzedni[k];
		  }
          
          wynik[j] += b[j]/A[j][j];
          
          if( kryterium & XN_MIN_XPOP && max( fabs( wynik[0] - poprzedni[0]), fabs( wynik[1] - poprzedni[1]) ) < tolx )
           break;
		  poprzedni[j] = wynik[j];
    }
    
	printf("Iteracja %d\t x1 = %g\t x2 = %g\t", iter, wynik[0], wynik[1] );
    printf("Residuum: %g\n", residuum( wynik ) );

    if ( kryterium & MAX_IT && iter >= liczba_iteracji )
     break;
    if( kryterium & RESIDUUM && fabs(residuum(wynik)) < tolf )
     break;
	}
}

//------------------ Metoda SOR -------------------
void ukladRownan::metodaSOR(double parametr)
{
 double wynik[2] = {przyblizenie_poczatkowe[0], przyblizenie_poczatkowe[1]};
 double poprzedni[2] = {przyblizenie_poczatkowe[0], przyblizenie_poczatkowe[1]};
 int iter = 0;
 double suma1 = 0 ,suma2 = 0;
 
 cout<<"Metoda SOR\n"<<endl;
 
	//wyliczanie przyblizen rozwiazan
	while(1)
	{
     ++iter;
		for(int j=0; j<n; j++)
		{
          wynik[j] = 0;
          suma1 = 0; suma2 = 0;
          
          for(int k=0; k<n; k++)
		  {
            if(k<j)
             suma1 += A[j][k]*poprzedni[k];
            if(k>=j)
             suma2 += A[j][k]*poprzedni[k];
		  }
          
          wynik[j] = poprzedni[j] - ((parametr/A[j][j])*(suma1+suma2-b[j]));
          
          if( kryterium & XN_MIN_XPOP && max( fabs( wynik[0] - poprzedni[0]), fabs( wynik[1] - poprzedni[1]) ) < tolx )
           break;
           
          poprzedni[j] = wynik[j];
		 }
         
         printf("Iteracja %d\t x1 = %g\t x2 = %g\t", iter, wynik[0], wynik[1] );
         printf("Residuum: %g\n", residuum( wynik ) );

         if ( kryterium & MAX_IT && iter >= liczba_iteracji )
          break;
         if( kryterium & RESIDUUM && fabs(residuum(wynik)) < tolf )
          break;
	}
}

//------------------ Metoda Residium -------------------
double ukladRownan::residuum( double x[2] )
{
  double Ax[2];
  Ax[0] = A[0][0]*x[0] + A[0][1]*x[1];
  Ax[1] = A[1][0]*x[0] + A[1][1]*x[1];
  
  return max(fabs(Ax[0]-b[0]), fabs(Ax[1]-b[1]));
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
int main()
{
   ukladRownan x;                  //Tworzenie obiektu.
   
   //------------------ Wydrukowanie macierzy A i wektora b z zad. 1 -------------------
   x.drukuj();
   
   
   x.metodaJacobiego();
   //x.metodaGaussaSeidela();
   //x.metodaSOR();
   system("pause");
    
}

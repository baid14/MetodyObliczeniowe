#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;
double EPS=0.000001;
double **new_matrix(int n, int m)                        //tworzenie macierzy
{
	int i;
	double **a;
	a=new double*[n];
	for(i=0; i<n; i++) a[i]=new double[m];
	return a;
}
void delete_matrix(double **a, int n)                     //usuwanie macierzy
{
	int i;
	for(i=n-1; i>=0; i--) delete []a[i];
	delete []a;
}
void wypelnij_macierz(double **a, int n)                               //wypelnianie macierzy
{
	int i,j;
	ifstream dane;
	dane.open("7.txt");
	for(i=0; i<n; i++)
	{
		for(j=0; j<n; j++)
		{
			dane>>a[i][j];
		}
	}
	dane.close();
}
void wyswietl_macierz(double **a, int n)                  //wyswietlanie macierzy
{
	int i,j;
	for(i=0; i<n; i++)
	{
		for(j=0; j<n; j++)
		{	
			cout<<a[i][j]<<" \t" ;
		}
		cout<<endl;
	}
}
void wyswietl_wektor(double *b,int n)            //wyswietlanie wektora 
{
     int i;
     for(i=0; i<n; i++)
     {
              cout<<b[i]<<endl;
     }
}
double norma(double *p)                         //norma max z wektora
{
	if(fabs(p[0])<fabs(p[1])) return fabs(p[1]);
	else return fabs(p[0]);
}
double* res(double **a, double *x, double *b, int n) //obliczanie wartsci wektora residualnego
{
     double suma=0;
     double *c=new double[n];
     int i,j;
     for(i=0; i<n; i++)
     {
              for(j=0; j<n; j++)
              {
                       suma=suma+a[i][j]*x[j];
              }
              c[i]=suma-b[i];
              suma=0;
     }
     return c;
} 
void jacobi(double **macierz, double *b, double *x,  int n)  //funkcja stosujaca metode iteracyjna Jacobiego
{
	int i, j, l=0;
	double *x1=new double[n];                          //wektory do przechowywania poprzednich przyblizen rozwiazan
	double *xnew;
	double *xold;
	double t=0;                                        //zmienna pomocnicza do obliczania sumy elementow macierzy i wektora b
 	while(1)
	{
            
                if(l==(l/2)*2)                           //warunek odpowiadajacy za zamiane nowych i poprzednich rozwiazan
                {
                               xnew=x1; 
                               xold=x;
                }
		        else 
                {
                     xnew=x; 
                     xold=x1;
                }
            for(i=0; i<n; i++)                          
            {
			for(j=0; j<n; j++)
			{
				if(j!=i) t=t+macierz[i][j]*xold[j];         //sumujemy elementy nie znajdujace sie na przekatnej 
			}
			xnew[i]=(1.0/macierz[i][i])*(b[i]-t);          //nowa wartosc pierwiastka jako rozw. rownania z macierza diagonalna
			t=0;
		}
		if(fabs((norma(res(macierz,xnew,b,2))))<EPS) break;  //warunek zakonczenia iteracji wartosci wektora residualnego do zera
		if(fabs((norma(xnew)-norma(xold)))<EPS) break;      //warunek zakonczenia iteracji normy przyblizen pierwiastkow do zera
		if(l>53) break;                                     //warunek zakonczenia iteracji skonczona ilosc iteracji
		l++;
		cout<<l<<" iteracja:"<<endl;
		wyswietl_wektor(xnew,2);
        cout<<"Wartosc normy z wektora residualnego                 wartosc normy z przyblizen pierwiastkow"<<endl;
        cout<<fabs(norma(res(macierz,xnew,b,2)));                  
        cout<<"                                             "<<fabs((norma(xnew)-norma(xold)))<<endl<<endl;                         
    }
}
void seidel(double **macierz, double *b, double *x, int n)  //funkcja stosujaca metoda Gaussa-Seidela
{
	int i,j, l=1;
	double *c=new double[n];                                            //wektor do przechowywania poprzednich przyblizen x;
	double t=0;                                                        //zmienna pomocnicza
	while(1)
	{
		for(i=0; i<n; i++)
		{
			c[i]=x[i];                                           //przechowanie poprzedniej wartosci pierwiastka
			for(j=0; j<=i-1; j++)
			{
				t=t+macierz[i][j]*x[j];                //suma elementow dolnej macierzy przy rozwiazywaniu ukladu z dolna macierza
			}	
			for(j=i+1; j<n; j++)
			{
				t=t+macierz[i][j]*x[j];                //suma elementow gornej macierzy
			}
			x[i]=(1.0/macierz[i][i])*(b[i]-t);        //wyznaczenie nowego przyblizenia x korzystajac z ukladu z dolna macierza
			t=0;
		}
	if(fabs((norma(res(macierz,x,b,2))))<EPS) break;  //warunek zakonczenia iteracji wartosci wektora residualnego do zera
	if(fabs((norma(x)-norma(c)))<EPS) break;      //warunek zakonczenia iteracji normy przyblizen pierwiastkow do zera
	if(l>53) break;                                     //warunek zakonczenia iteracji skonczona ilosc iteracji
	cout<<l<<" iteracja:"<<endl;
	wyswietl_wektor(x,2);
    cout<<"Wartosc normy z wektora residualnego                 wartosc normy z przyblizen pierwiastkow"<<endl;
    cout<<fabs(norma(res(macierz,x,b,2)));                  
    cout<<"                                             "<<fabs((norma(x)-norma(c)))<<endl<<endl;     
    l++;
	}
}
void sor(double **macierz, double *b, double *x, double omega, int n) //funkcja stosujaca metode SOR
{
    int i,j, l=1;
	double *c=new double[n];                                //wektor do przechowywania poprzednich wartosci pierwiastka
	double t=0;                                             //zmienna pomocnicza
	while(1)
	{
		for(i=0; i<n; i++)
		{
			c[i]=x[i];                                     //przechowanie poprzedniej wartosci pierwiastka
			for(j=0; j<=i-1; j++)
			{
				t=t+macierz[i][j]*x[j];                   //suma elementow dolnej macierzy pomnozona przez wektor 
			}	
			for(j=i; j<n; j++)
			{
                       if(j==i) t=t+(1-(1.0/omega))*macierz[i][i]*x[j];    //elementy przekatnej pomnozone przez parametr
                       else t=t+macierz[i][j]*x[j];                       //suma elementow gornej macierzy pomnoz przez wektor
			}
			x[i]=(1.0/macierz[i][i])*(b[i]-t)*omega;                  //rozwiazanie ukladu z dolna macierza
			t=0;
		}
	if(fabs((norma(res(macierz,x,b,2))))<EPS) break;  //warunek zakonczenia iteracji wartosci wektora residualnego do zera
	if(fabs((norma(x)-norma(c)))<EPS) break;      //warunek zakonczenia iteracji normy przyblizen pierwiastkow do zera
	if(l>53) break;                                     //warunek zakonczenia iteracji skonczona ilosc iteracji
	cout<<l<<" iteracja:"<<endl;
	wyswietl_wektor(x,2);
    cout<<"Wartosc normy z wektora residualnego                 wartosc normy z przyblizen pierwiastkow"<<endl;
    cout<<fabs(norma(res(macierz,x,b,2)));                  
    cout<<"                                             "<<fabs((norma(x)-norma(c)))<<endl<<endl;     
    l++;
	}
}
main()
{
      double **macierz;
      macierz=new_matrix(2,2);
      double b[2]={5.0, 6.0};
      double x[2]={1.0, 1.0};
      int a;
      wypelnij_macierz(macierz, 2);
      cout<<"Ktora metode chcesz zastosowac?"<<endl;
      cout<<"1. Jacobiego"<<endl<<"2. Gaussa-Seidela"<<endl<<"3. SOR"<<endl;
      cin>>a;
      switch(a)
      {
               case 1:
                    jacobi(macierz,b,x,2);
                    break;
               case 2:
                    seidel(macierz,b,x,2);
                    break;
               case 3:
                    sor(macierz,b,x,0.5,2);
                    break;
      }
      system("PAUSE");
      return 0;
}

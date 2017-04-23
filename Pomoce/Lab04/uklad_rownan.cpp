#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>

double rownanie_1(double x, double y)
{
       return x*x+y*y-0.5;
}
double rownanie_2(double x, double y)
{
       return 2*x-y*y*y;
}
double rownanie_1_ppo_x(double x)
{
       return 2*x;
}
double rownanie_1_ppo_y(double y)
{
       return 2*y;
}
double rownanie_2_ppo_x(double x)
{
       return 2;
}
double rownanie_2_ppo_y(double y)
{
       return -3*y*y;
}

using namespace std;

int main()
{
  int i=0;
  double   x[2],x1[2],wart_f[2],macierz_poch[2][2],d[2],wyznacznik,eps;

  cout<<"Podaj dokladnosc: ";
  cin>>eps;
  cout<<"Podaj punkt startowy (x,y): ";
  cin>>x[0]>>x[1];
  
  x1[0]=x[0]+1;
  x1[1]=x[1]+1;
  
  cout<<"\nObliczanie pierwiastka ukladu rownan - metoda Newtona\nx^2+y^2=0.5, 2*x-y^3=0\n";

  for(i=0;i<200 && fabs(x1[0]-x[0])>eps && fabs(x1[1]-x[1])>eps ;i++)
  {
	  x1[0]=x[0];
	  x1[1]=x[1];
  
	  //macierz pochodnych
	  macierz_poch[0][0]=rownanie_1_ppo_x(x[0]); macierz_poch[0][1]=rownanie_1_ppo_y(x[1]);
	  macierz_poch[1][0]=rownanie_2_ppo_x(x[0]); macierz_poch[1][1]=rownanie_2_ppo_y(x[1]);
  
	  //wektor f1 (wartosc rownania pierwszego i drugiego)
	  wart_f[0]=rownanie_1(x[0],x[1]);
	  wart_f[1]=rownanie_2(x[0],x[1]);
  
	  //wektor d
	  wyznacznik = macierz_poch[1][0] * macierz_poch[0][1] - macierz_poch[1][1] * macierz_poch[0][0];
	  if(fabs(wyznacznik) < 0.00000000000001)
		cout << "Brak rozwiazania";
	  else
	  {
		//poprawka ze wzgledu na X
		d[0] = (macierz_poch[1][1] * (-wart_f[0]) - (-wart_f[1]) * macierz_poch[0][1]) / wyznacznik;
		//poprawka ze wzgledu na Y
		d[1] = ((-wart_f[1]) * macierz_poch[0][0] - macierz_poch[1][0] * (-wart_f[0])) / wyznacznik;
	  }
	  printf("\niteracja: %d\tx = %.10lf\ty = %.10lf",i,x1[0],x1[1]);
  
	  //nowe przyblizenie
	  x[0]=x1[0]-d[0];
	  x[1]=x1[1]-d[1];
  }
  
  printf("\n\nWyniki:\nx = %.10lf\ny = %.10lf\n",x1[0],x1[1]);
  //cout<<"\nWyniki:\nx = "<<x[0]<<"\ny = "<<x[1]<<endl<<endl;
  cout<<"iteracji: "<<i<<endl<<endl;

  return 0;
} 

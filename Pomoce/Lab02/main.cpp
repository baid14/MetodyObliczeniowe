#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>

double r1(double x, double y)
{
       return x*x+y*y-8;
}
double r2(double x, double y)
{
       return y*y*y+4*x;
}
double r1pox(double x)
{
       return 2*x;
}
double r1poy(double y)
{
       return 2*y;
}
double r2pox(double x)
{
       return 4;
}
double r2poy(double y)
{
       return 3*y*y;
}

/*double r1(double x, double y)
{
       return x*x-y+1;
}
double r2(double x, double y)
{
       return 3*cos(x)-y;
}
double r1pox(double x)
{
       return 2*x;
}
double r1poy(double y)
{
       return -1;
}
double r2pox(double x)
{
       return -3*sin(x);
}
double r2poy(double y)
{
       return -1;
}*/


using namespace std;

int main()
{
  int i=0;
  double   x[2],x1[2],f1[2],mp[2][2],d[2],m,eps;

  cout<<"Podaj dokladnosc: ";
  cin>>eps;
  cout<<"Podaj punkt startowy (x,y): ";
  cin>>x[0]>>x[1];
  
  x1[0]=x[0]+1;
  x1[1]=x[1]+1;
  
  cout<<"--------------------------------";
  for(i=0;i<200 && fabs(x1[0]-x[0])>eps && fabs(x1[1]-x[1])>eps ;i++)
  {
  x1[0]=x[0];
  x1[1]=x[1];
  
  //macierz pochodnych
  mp[0][0]=r1pox(x[0]); mp[0][1]=r1poy(x[1]);
  mp[1][0]=r2pox(x[0]); mp[1][1]=r2poy(x[1]);
  
  //wektor f1
  f1[0]=r1(x[0],x[1]);
  f1[1]=r2(x[0],x[1]);
  
  //wektor d
  //a = mp[0][0]
  //b = mp[0][1]
  //c = f1[0]
  //d = mp[1][0]
  //e = mp[1][1]
  //f = f1[1]
  m = mp[1][0] * mp[0][1] - mp[1][1] * mp[0][0];
  if(fabs(m) < 0.00000000000001)
    cout << "Brak rozwiazania";
  else
  {
    d[0] = (mp[1][1] * (-f1[0]) - (-f1[1]) * mp[0][1]) / m;
    d[1] = ((-f1[1]) * mp[0][0] - mp[1][0] * (-f1[0])) / m;
  }
  printf("\niteracja: %d\tx = %.10lf\ty = %.10lf",i,x1[0],x1[1]);
  
  //nowe przyblizenie
  x[0]=x1[0]-d[0];
  x[1]=x1[1]-d[1];
  }
  
  printf("\n-------------------------\nWyniki:\t\t\t|\nx = %.10lf\t|\ny = %.10lf\t|\n\t\t\t|\n",x1[0],x1[1]);
  //cout<<"\nWyniki:\nx = "<<x[0]<<"\ny = "<<x[1]<<endl<<endl;
  cout<<"-------------------------\niteracji: "<<i<<endl<<endl;

  system("pause");

  return 0;
} 

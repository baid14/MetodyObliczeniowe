#include <iostream>
#include <cmath>
using namespace std;
double roznica=0;
double h=0.1;           //poczatkowa wartosc kroku h
void forward(double x0, double x1, double x2)    //funkcja stosujaca roznice wprzod dla pkt poczatkowych
{
     roznica=(-1.5*exp(x0)+2*exp(x1)-0.5*exp(x2))/(h);
}
void backward(double xn2, double xn1, double xn)     //funkcja stosujaca roznice wstecz dla pkt koncowych
{
     roznica=(0.5*exp(xn2)-2*exp(xn1)+1.5*exp(xn))/(h);
}
void central(double xi1, double xi2)       //funkcja stosujaca roznice centralna
{
     roznica=(exp(xi2)-exp(xi1))/(2*h);
}
main()
{
      double blad;                       //wartosc przechowujaca blad danej roznicy
      double dokladnosc;                //wartosc przechowujaca dokladnosc danej roznicy
      double stop;                      //zmienna potrzebna do liczenia epsylonu maszynowego
      while(1)
      {
      cout<<"h"<<"           | "<<"poczatek   |  blad      | dokladnosc |"<<"koniec      | blad       |dokladnosc  |"<<"centralna   |";    
      cout<<"blad        "<<"| dokladnosc"<<endl;
      cout.setf(ios_base::scientific);          //wyswietlanie wynikow w postaci naukowej
      forward(0,h,2*h);
      cout<<h<<"|"<<roznica;
      blad=fabs((roznica-exp(0))/(exp(0)));    //liczenie bledu 
      dokladnosc=fabs(roznica-exp(0));         //liczenie dokladnosci
      cout<<"|"<<blad<<"|"<<dokladnosc;
      backward(1-2*h,1-h,1);
      cout<<"|"<<roznica;
      blad=fabs((roznica-exp(1))/(exp(1)));
      dokladnosc=fabs(roznica-exp(1));
      cout<<"|"<<blad<<"|"<<dokladnosc;
      central(0.5-h,0.5+h);
      cout<<"|"<<roznica;
      blad=fabs((roznica-exp(0.5))/(exp(0.5)));
      dokladnosc=fabs(roznica-exp(0.5));
      cout<<"|"<<blad<<"|"<<dokladnosc<<endl;
cout<<"------------|------------|------------|------------|------------|------------|------------|------------|------------|------------"<<endl;
      stop=h+1.0;     
      if(stop<=1.0)break;               //osiagniecie wartosci h rownej epsylonowi maszynemowemu to zatrzymanie obliczen
      h=h*0.1;                          //nastepny krok h w skali logarytmicznej
      }
      system("PAUSE");
      return 0;
}

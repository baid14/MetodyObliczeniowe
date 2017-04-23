#include <stdio.h>
#include<windows.h>
#include<math.h>

/*

Metody Obliczeniowe
By Wahuu

*/


////////////////////////////////  Zmienne globalne  //////////////////////////////////

int max_il_krokow = 50;
////////////////////////////////  Dwie funkcje  //////////////////////////////////

double rownfun(double x0)
{
   double w = cos(x0/2.0) +1.0;   //Rownanie funkcji 1;
   //double w = 1.0/tan(x0);     //Rownanie funkcji 2;
   return w;
}

double fun(double x)
{
   double w = x-cos(x/2)-1;         //Funkcja 1
   //double w = x-(1.0/tan(x));    //Funkcja 2
   return w;
}

double dfun(double x)  //pochodne funkcji
{
   double w = 1+(0.5*sin(x/2));           //Funkcja 1 - pochodna
   //double w = 1+(1/(sin(x)*sin(x)));   //Funkcja 2 - pochodna
   return w;
}


/////////////  Deklaracje funkcji + main + deklaracje globalnych  ////////////////

double EPS=10e-10, fEPS=10e-10;  //dokladnosc, dokladnosc wartosci
void picard();
void bisekcja();
void newton();
void sieczne();

int main(void){
	
	printf("METODY OBLICZENIOWE Zad 3\n\n");
	printf("****************************************************\n\n");
	printf("Nacisnij dowolny klawisz aby kontynuowac\n\n");
	getch();
	
	picard();
	bisekcja();
	newton();
	sieczne();
	
	return 0;
}


//////////////////////////////////  Picard  //////////////////////////////////////

void picard()
{
   printf("\n-----------------------  Picard  -----------------------\n\n");
   int i=0;//zmienna iteracyjna
   double  x=0.5;//wartosc poczatkowa
   double x0=x+0.1;//Zapamietanie wartosci poczatkowej
   
   while(i++<=max_il_krokow  &&  fabs(x-x0)>EPS  &&  fabs(fun(x))>fEPS)
   {
      x0 = x;
      x = rownfun(x);
      printf("i = %3d \t %lf\n", i, x);
   }
   printf("-------------------  Picard - koniec  --------------------\n");
   getch();
}

//////////////////////////////////  Bisekcja  ////////////////////////////////////

void bisekcja()
{
   printf("\n-----------------------  Bisekcja  -----------------------\n\n");
   double a=0.5, b=2.0,c;
   int i=0;
   while(i++<=max_il_krokow  &&  fun(a)!=0  &&  fun(b)!=0  &&  fabs((b-a)/2.)>=EPS  &&  fabs(fun(c))>fEPS)
   {
      c=(a+b)/2;
      if((fun(a)<0 && fun(c)>0) || (fun(a)>0 && fun(c)<0))
      {
         b=c;
      }
      else
      {
         if((fun(c)<0 && fun(b)>0) || (fun(c)>0 && fun(b)<0))
         {
            a=c;
         }
         else
         {
            printf("Brak wyniku"); break;
         }
      }
      printf("i = %3d \t %lf\n", i, c);
   }
   printf("------------------  Bisekcja - koniec  ------------------\n");
   getch();
}

/////////////////////////////////  Newton  //////////////////////////////////////


void newton()
{
    printf("\n-----------------------  Newton  ------------------------\n\n");
   int i=0;//zmienna iteracyjna
   double  x=0.5;//wart poczatkowa
   double x0=x+0.1;//Zapamietanie wartosci poczatkowej
   
   while(i++<=max_il_krokow  &&  fabs(x-x0)>EPS  &&  fabs(fun(x))>fEPS)
   {
      x0 = x;
      x = x - (fun(x)/dfun(x));
      printf("i = %3d \t %lf\n", i, x);
   }
   printf("-------------------  Newton - koniec  -------------------\n");
   getch();
}

////////////////////////////////  Sieczne  ///////////////////////////////////////

void sieczne()
{
   printf("\n-----------------------  Sieczne  -----------------------\n\n");
   int i=0;//zmienna iteracyjna
   double x1=0.4;//wart poczatkowa
   double x3=0.6, x2 = 0.5;

   while(i++<=max_il_krokow  &&  fabs(x3-x2)>EPS  &&  fabs(fun(x3))>fEPS)
   {
      x1 = x2;
      x2 = x3;
      x3 = x2 - (fun(x2)*(x2-x1))/(fun(x2)-fun(x1));
      printf("i = %3d \t %lf\n", i, x3);
   }
   printf("------------------  Sieczne - koniec  -------------------\n");
   getch();
}

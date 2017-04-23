#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXITER 200

double f(double x)
{
       return x-atan(x);
       //return x-sin(x/3);
       //return x-tan(x/2);
}

double ff(double x)
{
       return 1-1/(1+x*x);
       //return 1-1/3*cos(x/3);
       //return 1-0.5/(cos(x/2)*cos(x/2));
}

double fe(double x)
{
       return atan(x);
       //return sin(x/3);
       //return tan(x/2);
}

double ffe(double x)
{
       return 1/(1+x*x);
       //return 1/3*cos(x/3);
       //return 0.5/(cos(x/2)*cos(x/2));
}

int main(int argc, char *argv[])
{
    int i,ii,menu,exit=0,test=0;
    double a,b,c,c2,d,x0,f0,eps;
    
    while(exit==0){
    printf("\nWybierz metode: \n 1 -metoda bisekcji\n");
    printf(" 2 -metoda iteracji prostych\n 3 -metoda Regula falsi\n");
    printf(" 4 -metoda siecznych\n 5 -metoda Newtona\n 6 -wyjscie\n\n");
    scanf("%d",&menu);
    
    switch(menu)
    {          
          case 1:
               printf("\nMetoda bisekcji\n");
               printf("\nPodaj dokladnosc: ");
               scanf("%lf",&eps);
               printf("Podaj przedzial <a,b>: ");
               scanf("%lf %lf",&a,&b);
               if (f(a)*f(b) > 0)
               {
                   printf("Niewlasciwy przedzial!\n");
                   break;
               }
               ii=0;
               for(i=0;i<MAXITER && (b-a)/2>eps;i++)
               {
                    c = (a+b)/2;
                    if (f(c)==0) break;
                    if (f(a)*f(c)<0) b = c; // pierwiastek jest w podprzedziale (a;c)
                    else a = c; // pierwiastek jest w podprzedziale (c;b)
                    //printf("\niteracja: %d\tx = %.10lf",i,c);
               }
               printf("\n\nx = %.10lf +/- %.10lf\niteracji: %d\n",c,(b-a)/2,i);
               break;
          case 2:
               printf("\nMetoda iteracji prostych\n");
               printf("\nPodaj dokladnosc: ");
               scanf("%lf",&eps);
               printf("Podaj punkt startowy: ");
               scanf("%lf",&a);
               b=a+1;
               if(ffe(0)>=1) 
               {
                                 printf("Metoda rozbiezna dla tej funkcji\n");
                                 break;
               }
               for(i=0;i<MAXITER && fabs(b - a)>eps;i++)
               {
                    b=a;
                    a=fe(a);
                    //printf("\niteracja: %d\tx = %.10lf",i,a);
               }
               if(!test) printf("\n\nx = %.10lf\niteracji: %d\n",a,i);
               break;
          case 3:
               printf("\nMetoda Regula falsi\n");
               printf("\nPodaj dokladnosc: ");
               scanf("%lf",&eps);
               printf("Podaj przedzial <a,b>: ");
               scanf("%lf %lf",&a,&b);
               if (f(a)*f(b) > 0)
               {
                   printf("Niewlasciwy przedzial!\n");
                   break;
               }
               c = (a*f(b)-b*f(a))/(f(b)-f(a));
               c2+c+1;
               for(i=0;i<MAXITER && fabs(c-c2)>eps;i++)
               {
                    c2=c;
                    if (f(c)==0) break;
                    if (f(a)*f(c)<0) b = c; // pierwiastek jest w podprzedziale (a;c)
                    else a = c; // pierwiastek jest w podprzedziale (c;b)
                    c = (a*f(b)-b*f(a))/(f(b)-f(a));
                    //printf("\niteracja: %d\tx = %.10lf",i,c);
               }
               printf("\n\nx = %.10lf\niteracji: %d\n",c,i);
               break;
          case 4:
               printf("\nMetoda siecznych\n");
               printf("\nPodaj dokladnosc: ");
               scanf("%lf",&eps);
               printf("Podaj dwa wstepne punkty x1 i x2: ");
               scanf("%lf %lf",&a,&b);
               if(a==b)
               {
                    printf("Punkty nie moga byc rowne!\n");
                    break;
               }
               c = f(a); d = f(b); i = 1;
               for(i=0;i<MAXITER;i++)
               {
                    x0 = a - c * (a - b) / (c - d);
                    f0 = f(x0);
                    b = a; d = c;
                    a = x0; c = f0;
                    if(fabs(a - b) < eps/1000) break;
                    //printf("\niteracja: %d\tx = %.10lf",i,x0);
               }
               if(i) printf("\n\nx = %.10lf\niteracji: %d\n",x0,i);
               break;
          case 5:
               printf("\nMetoda Newtona\n");
               printf("\nPodaj dokladnosc: ");
               scanf("%lf",&eps);
               printf("Podaj punkt startowy: ");
               scanf("%lf",&a);
               b = a + 1;
               for(i=0;i<MAXITER;i++)
               {
                     b = a;
                     a = b - f(b) / ff(b);
                     if(fabs(a - b) < eps/1000) break;
                     //printf("\niteracja: %d\tx = %.10lf",i,a);
               }
               if(i) printf("\n\nx = %.10lf\niteracji: %d\n",a,i);
               break;
          case 6:
               exit=1;
               break;
          default:
               printf("\nBrak opcji!\n");
               break;
    }
    }
  
  //system("PAUSE");	
  return 0;
}

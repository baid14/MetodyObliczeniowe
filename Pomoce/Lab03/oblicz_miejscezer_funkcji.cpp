#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_ITERACJI 200
#define ROWNANIE 1

double f(double x)
{
#if ROWNANIE == 1
       //return x-atan(x);
       return pow( cos( x / 4 ), 2 ) - x;
#elif ROWNANIE == 2
       return x-sin(x/3);
#else
       return x-tan(x/2);
#endif
}

double f_pochodna(double x)
{
#if ROWNANIE == 1
       //return 1-1/(1+x*x);
       return -0.25 * sin( x / 4.0 );
#elif ROWNANIE == 2
       return 1-1/3*cos(x/3);
#else
       return 1-0.5/(cos(x/2)*cos(x/2));
#endif
}

double f_picarda(double x)
{
#if ROWNANIE == 1
       return atan(x);
#elif ROWNANIE == 2
       return sin(x/3);
#else
       return tan(x/2);
#endif
}

double f_picarda_pochodna(double x)
{
#if ROWNANIE == 1
       return 1/(1+x*x);
#elif ROWNANIE == 2
       return 1/3*cos(x/3);
#else
       return 0.5/(cos(x/2)*cos(x/2));
#endif
}

int main()
{
	int wybor, koniec = 0, i;
	double eps, a, b,c,c2;

	while(koniec==0)
	{
		printf("\nMENU\nWybierz metode ktra chcesz zastosowac: \n 1 - Metoda bisekcji\n");
		printf(" 2 - Metoda iteracji prostych (metoda Picarda)\n 3 - Regula falsi\n");
		printf(" 4 - Metoda siecznych\n 5 - Metoda Newtona\n 6 - Wyjscie\n\n");
		scanf("%d",&wybor);

		switch(wybor)
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

               for(i=0;i<MAX_ITERACJI && (b-a)/2>eps;i++)
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

               if(f_picarda_pochodna(0)>=1)
               {
                                 printf("Metoda rozbiezna dla tej funkcji\n");
                                 break;
               }

               for(i=0;i<MAX_ITERACJI && fabs(b - a)>eps;i++)
               {
                    b=a;
                    a=f_picarda(a);
                    //printf("\niteracja: %d\tx = %.10lf",i,a);
               }

               printf("\n\nx = %.10lf\niteracji: %d\n",a,i);
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
               c2=c+1;

               for(i=0;i<MAX_ITERACJI && fabs(c-c2)>eps;i++)
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

			   double d, x0, f0;
               c = f(a); d = f(b); i = 1;
               for(i=0;i<MAX_ITERACJI;i++)
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

               for(i=0;i<MAX_ITERACJI;i++)
			   {
                     b = a;
                     a = b - f(b) / f_pochodna(b);
                     if(fabs(a - b) < eps/1000) break;
                     //printf("\niteracja: %d\tx = %.10lf",i,a);
               }

               if(i) printf("\n\nx = %.10lf\niteracji: %d\n",a,i);
               break;

          case 6:
               koniec=1;
               break;

          default:
               printf("\nBrak opcji!\n");
               break;
		}
	}
}

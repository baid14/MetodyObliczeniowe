#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define DEBUG 0
const double EPS = 0.0000000001;
const int MAX = 100;
int i;
//Funkcja 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double funkcja1(double x)
{
    return 1-sin(x/2)-x;
}

double funkcja1_pochodna(double x)
{
    return 1/2*cos(x/2)-1;
}

double funkcja1picard(double x)
{
    return 1-sin(x/2);
}

double funkcja1picard_pochodna(double x)
{
    return 1/2*cos(x/2);
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


// Funkcja 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double funkcja2(double x)
{
    return 1-tan(x)-x;
}

double funkcja2_pochodna(double x)
{
    return 1/(cos(x)*cos(x))-1;
}

double funkcja2picard(double x)
{
    return 1-tan(x);
}

double funkcja2picard_pochodna(double x)
{
    return 1/(cos(x)*cos(x));
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Metoda Bisekcji %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double bisekcja(double (*f)(double), double a, double b)
{
    double c=0;

    printf("Metoda bisekcji\n");
    if(f(a)*f(b) < 0.0)
        {

            for(i=0; i<MAX && fabs(b-a)/2>EPS; i++)
            {
                c = (a+b)/2.0;
                if (f(c)==0.0)
                    break;
                if((f(a) * f(c)) < 0)
                    b = c; // Pierwiastek znajduje sie w podprzedziale (a;c)
                    else
                        a = c; // Pierwiastek znajduje sie w podprzedziale (c;b)
                if(DEBUG) printf("Iteracja[%d]: \tx = %.20lf\n",i+1,c);
            }
        }

return c;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//Metoda Newtona %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double newton(double (*f)(double), double (*f_pochodna)(double),double a)
{
	double b;

    printf("Metoda Newtona\n");
	for(i=0; i<MAX; i++)
	{
	   b = a;
	   a = b - f(b)/f_pochodna(b);
	   if(fabs(a-b)<EPS)
            break;
	   if(DEBUG) printf("Iteracja[%d]: \tx = %.20lf\n",i+1,a);
	}
return a;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//Metoda Picarda %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double picard(double (*f)(double), double (*f_poch)(double), double a)
{
    if(f_poch(0)>=1.0)
    {
		printf("Metoda rozbiezna dla tej funkcji\n");
    }
	else
	{
    double b;
    b=a+1.0;
	for(i=0; i<MAX && fabs((a-b)/2)>EPS; i++)
		{
			//b jest poprzedni¹ wartoœci¹ a
			b=a;
			// liczymy wartoœæ funkcji dla a ( na prostej y=x, wiêc a, która jest wartoœci¹ funkcji traktujemy jako argument)
			a=f(a);

			//warunek zatrzymania
			if(!(fabs(b - a)>0))
				break;
            if(DEBUG) printf("Iteracja[%d]: \tx = %.20lf\n",i+1,a);

		}
	}
    return a;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//Metoda siecznych %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
double sieczne(double (*f)(double), double a, double b)
{
	double c, d, x0, f0;
	if(a==b)
    {
		printf("Punkty nie moga byc rowne!\n");
    }
	else
	{
		c = f(a);
		d = f(b);

		for(i=0;i<MAX && fabs(a - b)> EPS;i++)
		{
			x0 = a - c * (a - b) / (c - d);
			f0 = f(x0);
			b = a; d = c;
			a = x0; c = f0;
			if(fabs(a - b) < EPS) break;
			if(DEBUG) printf("Iteracja[%d]: \tx = %.20lf\n",i+1,x0);
		 }
	}

    return x0;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


int main()
{
    printf("\n\tMETODA BISEKCJI\nFUNKCJA 1\nPierwiastek wynosi: %.20lf \nz dokladnoscia %.e\n-----------------------\n",bisekcja(funkcja1,-2,2), EPS);
    printf("\n\tMETODA BISEKCJI\nFUNKCJA 2\nPierwiastek wynosi: %.20lf \nz dokladnoscia %.e\n-----------------------\n",bisekcja(funkcja2,-5,6), EPS);
    printf("\n\tMETODA NEWTONA\nFUNKCJA 1\nPierwiastek wynosi: %.20lf \nz dokladnoscia %.e\n-----------------------\n",newton(funkcja1,funkcja1_pochodna,10), EPS);
    printf("\n\tMETODA NEWTONA\nFUNKCJA 2\nPierwiastek wynosi: %.20lf \nz dokladnoscia %.e\n-----------------------\n",newton(funkcja2,funkcja2_pochodna,1), EPS);
    printf("\n\tMETODA PICARDA\nFUNKCJA 1\nPierwiastek wynosi: %.20lf \nz dokladnoscia %.e\n-----------------------\n",picard(funkcja1picard,funkcja1picard_pochodna,1), EPS);
    printf("\n\tMETODA PICARDA\nFUNKCJA 2\nPierwiastek wynosi: %.20lf \nz dokladnoscia %.e\n-----------------------\n",picard(funkcja2picard,funkcja2picard_pochodna,2), EPS);
    printf("\n\tMETODA SIECZNYCH\nFUNKCJA 1\nPierwiastek wynosi: %.20lf \nz dokladnoscia %.e\n-----------------------\n",sieczne(funkcja1,1,5), EPS);
    printf("\n\tMETODA SIECZNYCH\nFUNKCJA 2\nPierwiastek wynosi: %.20lf \nz dokladnoscia %.e\n-----------------------\n",sieczne(funkcja2,2,10), EPS);
    getch();
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ITERACJE 1000                   //Maksymalna ilosc iteracji
#define EPS 0.00000001                  //Dok³adnoœæ

enum
{
    PICARDA = 1,
    BISEKCJI = 2,
    NEWTONA = 3,
    SIECZNYCH = 4,
    FINISH = 5
};

//--------- Obliczanie wartoœci wybranej Funkcji ----------
double fun(double x, int n)
{
       if(n==1)
               return x+sin(x/2)-1;                  //Zad. 1.a
       if(n==2)
               return x+tan(x)-1;                    //Zad. 1.b
}

//--------- Obliczanie wartoœci pochodnej wybranej Funkcji ----------
double pochodna(double x, int n)
{
       if(n==1)
               return 1+cos(x/2)*0.5;                  //Zad. 1.a
       if(n==2)
               return 1+1/(cos(x)*cos(x));             //Zad. 1.b
}

//--------- Obliczanie wartoœci wybranej Funkcji w postaci x=f(x) ----------
double fun_picard(double x, int n)
{
       if(n==1)
               return 1-sin(x/2);                  //Zad. 1.a
       if(n==2)
               return 1-tan(x);                    //Zad. 1.b
}

//--------- Obliczanie wartoœci pochodnej wybranej Funkcji w postaci x=f(x) ----------
double pochodna_picard(double x, int n)
{
       if(n==1)
               return -cos(x/2)*0.5;                           //Zad. 1.a
       if(n==2)
               return -1.0/(cos(x)*cos(x));                    //Zad. 1.b
}

//--------- Obliczanie wartoœci pochodnej wybranej Funkcji w postaci do metody Newtona ----------
double do_newton(double x, int n)
{
       return x - fun(x, n)/pochodna(x, n);
}


/*****************************************************************************
******************************************************************************/

//--------- METODA PICARDA----------
double picard(double x, int n)
{
       int i=0;
       double xn_pop=x;
       double xn=fun_picard(x, n);
       
       while(i < ITERACJE && fabs(xn - xn_pop) > EPS )
       {
               xn_pop = xn;
               xn = fun_picard(xn, n);
               
               printf("*** Iteracja = %d, %.10lf\n",i, xn);
               ++i;
       }
       
       printf("Iteracja %d, %.10lf\n\n", i-1, xn);
}

//--------- METODA BISEKCJI----------
double bisekcja(int n, double a, double b)
{
       int i;
       double c;
       
       for( i=0; i < ITERACJE && (b-a)/2 > EPS; i++)
       {
            c=(a+b)/2;
            
            if( fun(c, n) == 0)
              break;
            
            if( fun(a, n) * fun(c, n) <0)
              b=c;                       //Pierwiastek znajduje sie w przedziale (a,c)
            else
              a=c;                       //Pierwiastek znajduje sie w przedziale (c,b)
              
            printf("***Iteracja = %d, %.10lf\n", i, c);
       }
       
       printf("Iteracja %d, %.10lf\n\n", i, c);
       
       return c;
}

//--------- METODA NEWTONA----------
double newton(double x, int n)
{
       int i=1;
       double xn1;
       double xn2 = do_newton(x, n);
       
       while (i < ITERACJE && fabs(xn1-xn2) >= EPS)
       {
             printf("*** Iteracja = %d, %.10lf\n", i, xn2);
             
             xn1 = xn2;
             xn2 = do_newton(xn2, n);
             ++i;
       }
       
       return xn2;       
}


//--------- METODA SIECZNYCH----------
double sieczne(int n, double a, double b)
{
       int i;
       double fa, fb, xn, f0;
       
       fa = fun(a, n);
       fb = fun(b, n);
       
       for(i=0; i < ITERACJE; i++)
       {
         xn = a-fa*(a-b)/(fa-fb);
         f0 = fun(xn, n);
         
         b = a;
         fb = fa;
         a = xn;
         fa = f0;
         
         if(fabs(a-b) < EPS)
           break;
           
         printf("***Iteracja = %d  x = %.10lf\n", i, xn);
       }
       
       if(i)
         printf("Iteracja: %d  x = %.10lf\n\n", i-1, xn);
}

/*****************************************************************************
******************************************************************************/

int main()
{
    int n;
    int nr_metody;
    double x;
    double a, b;
    
    while(nr_metody != 5)
    {
    printf("Ktora funkcje wybierasz (wprowadz numer):\n");
    printf("1: a) \n2: b)\nNumer funkcji: ");
    
    scanf("%d", &n);
    
    system("cls");
      printf("|-----------------------------------------------------------------------|\n");
      printf("|   ***                    MENU                                  ***    |\n");
      printf("|   ****************************************************************    |\n");
      printf("|   ****************************************************************    |\n");
      printf("|   *     Wybierz jedna z metod wybierajac jej numer               *    |\n");
      printf("|   ****************************************************************    |\n");
      printf("|   *      1 - metoda PICARDA                                      *    |\n");
      printf("|   *      2 - metoda BISEKCJI                                     *    |\n");
      printf("|   *      3 - metoda NEWTONA                                      *    |\n");
      printf("|   *      4 - metoda SIECZNYCH                                    *    |\n");
      printf("|   *      5 - Zakoncz dzialanie programu                          *    |\n");
      printf("|-----------------------------------------------------------------------|\n\n");
      
      printf("Numer: ");
      scanf("%d", &nr_metody);
      
      switch(nr_metody)
      {
        case PICARDA:
             printf("Podaj x: ");
             scanf("%lf",&x);
             
             picard(x, n); 
             break;
             
        case BISEKCJI:
             printf("Podaj krance przedzialu \n a = ");
             scanf("%lf", &a);
             printf("\n b = ");
             scanf("%lf", &b);
             
             if( fun(a, n) * fun(b, n) >0)
               {
                 printf("Podany przedzial jest niepoprawny!\n\n");
                 break;
               }
             
             bisekcja(n, a, b);
             break;
        
        case NEWTONA:
             printf("Podaj x: ");
             scanf("%lf", &x);
             
             newton(x, n);
             break;
             
        case SIECZNYCH:
             printf("Podaj dwa punkty startowe\n x1 = ");
             scanf("%lf", &a);
             printf("\n x2 = ");
             scanf("%lf", &b);
             
             sieczne(n, a, b);
             
             if(a == b)
             {
                  printf("Dwa punkty nie moga byc takie same!\n\n");
             }
             break;
             
        case FINISH:
                  exit(-10);
             
        default:
             printf("Podany numer jest nieprawidlowy!\n\n");
      }
      
    }
    
    system("pause");
      
}

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//dokladnosc 1
const double eps = 0.00001;
//dokladnosc 2
const double eps2 = 0.00001;
//siatka
const double dx = 0.00001;
//ilosc obrotow
const double loop = 52;
// x + sin(x/2) - 1 = 0
double f1(double x) {
    return (x + sin(x/2) - 1);}
// x = 1 - sin(x/2)
double f1_picard(double x) {
    return (1 - sin(x/2));}
// x + tan(x) - 1 = 0
double f2(double x) {
    return (x + tan(x) - 1);}
// x = 1 - tan(x)
double f2_picard(double x) {
    return (1 - tan(x));}
//metoda roznicy centralnej
double pochodna(double (*f)(double), double x, double dx) {
    return ((f)(x + dx) - f(x - dx)) / (2 * dx);}
/*
Metoda wyznaczajaca kolejne przyblizenia przy pomocy siecznych w danym punkcie.
*/
double sieczne(double (*f)(double x), double x1, double x2) {
    double f1 = f(x1), f2 = f(x2), x0, f0;
    int i = 1;
    printf("Metoda siecznych\n");
    do {
        x0 = x1 - f1 * (x1 - x2) / (f1 - f2);
        f0 = f(x0);
        x2 = x1;
        f2 = f1;
        x1 = x0;
        f1 = f0;
        printf("%d x= %lf\n",i,x0);
        i++;
    }
    while(i <= loop && (fabs(x1 - x2) > eps) && fabs(f1 - f2) > eps && fabs(f0) >= eps2);
    return x0;
}
/*Metoda wyznaczajaca kolejne przyblizenia przy pomocy metody picarda.
Rzutowanie wartosci funkcji w punkcie na prosta y = x.
Wartosc funkcji w punkcie staje sie kolejnym przyblizeniem  pierwiastka.
*/
 double picard(double (f)(double), double x0) {   
    int i = 1;
    double f0 = f(x0);
    printf("Metoda Picarda\n");
    while(fabs(f0 - x0) > eps && i <= loop) {
        x0 = f0;
        f0 = f(x0);
        printf("%d x= %lf\n",i,x0);
        i++;
    }
    return x0;
}
/*
Metoda wyznaczajaca kolejne przyblizenia przy pomocy pochodnych w danym punkcie.
*/
 double newton(double (*f)(double), double x0) {
    double f1, x1 = x0 + 0.1, f0 = f(x0);
    int i = 1;
    printf("Metoda Newtona\n");
    do {
        f1 = pochodna(f, x0, dx);
        x1 = x0;
        x0 = x0 - f0 / f1;
        f0 = f(x0);
        printf("%d x= %lf\n",i,x0);
        i++;
    }
    while(i <= loop && fabs(x1 - x0) > eps && fabs(f1) >= eps2);
    return x0;
}
/*
Metoda bisekcji znajduje tylko jeden pierwiastek, stopnia nieparzystego.
*/
 double bisekcja(double(*f)(double), double a, double b) {
    int i = 1;
    double x0;
    printf("Metoda bisekcji\n");
    if (f(a) * f(b) < 0) {
        while(i < loop && fabs(b-a) >= eps && fabs(f(b) - f(a)) >= eps2) {
        x0 = (a + b) / 2.0;
        if (f(a) * f(x0) < 0 )
          b= x0;
        else
          a = x0;
          printf("%d x= %lf\n",i,x0);
        i++;
        }
        return x0;
    }
    else
        printf("Bledny przedzial.\n");
        return 0;
}
 
int main(void) {
    printf("(sin): x= %lf\n",picard(f1_picard,5));
    printf("(sin): x= %lf\n",bisekcja(f1,-5,5));
    printf("(sin): x= %lf\n",newton(f1,1));
    printf("(sin): x= %lf\n",sieczne(f1,1,5));
   
    printf("(tan): x= %lf\n",picard(f2_picard,5));
    printf("(tan): x= %lf\n",bisekcja(f2,-5,5));
    printf("(tan): x= %lf\n",newton(f2,1));
    printf("(tan): x= %lf\n",sieczne(f2,1,5));

    system("pause");
}

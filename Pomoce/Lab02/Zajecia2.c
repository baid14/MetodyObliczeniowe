#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main()
{
    double a, b, c;                 //Wspolczynniki rownania kwadratowego.
    double delta;
    double x0_szkolne, x1_szkolne;
    double x0_poprawne, x1_poprawne;
    double blad0_szkolny, blad1_szkolny;
    double blad_poprawne;

    double x0_znane=1e-8;              //Znane miejsca zerowe.
    double x1_znane=3.0;               //Znane miejsce zerowe.

    a=2.0;

    /*********** SPAWDZENIE, CZY ROWNANIE KWADRATOWE ***********/
    if(a == 0)
	{
        printf("To nie jest rownanie kwadratowe!\n");
        exit(1);
	}

	printf("-----> Znane miejsca zerowe: \n x0_znane = %lf \n x1_znane = %lf\n\n", x0_znane, x1_znane);

    //---------- WZORY VIETE'a ------------------------------
    /** Obliczanie wspolczynnikow 'b' oraz 'c' na podstawie wzorow Vietea **/

        /* -----> Suma pierwiastkow: x0 + x1 = -b/a  =>   b = -a(x0 + x1)  <---- */
        /* -----> Iloczyn pierwiastkow: x0 * x1 = c/a  =>   c = a * x0 * x1  <---- */

	b = -a * (x0_znane * x1_znane);
	c = a * x0_znane * x1_znane;


    printf("-----> Wspolczynniki rownania kwadratowego: \n %.101f \n%.101f \n%.101f\n\n", a, b, c);

    //---------- WZORY SZKOLNE ---------------------------------

    delta = (b * b) - (4 * a * c);

	x0_szkolne = (-b - sqrt(delta))/(2 * a);
	x1_szkolne = (-b + sqrt(delta))/(2 * a);

	printf("-----> Pierwiastki ze wzorow szkolnych: \n x0_szkolne = %.20lf\n x1_szkolne = %.20lf\n\n", x0_szkolne, x1_szkolne);


    //---------- BLEDY WZOROW SZKOLNYCH ---------------------------------

	blad0_szkolny = (x0_szkolne - x0_znane)/x0_znane;
	blad1_szkolny = (x1_szkolne - x1_znane)/x1_znane;


   //---------- WZOR POPRAWNY ---------------------------------

	x0_poprawne = c/(a * x1_szkolne);


   //---------- BLAD WZORU POPRAWNEGO ---------------------------------

	blad_poprawne = (x0_poprawne - x0_znane)/x0_znane;

    system("pause");
}

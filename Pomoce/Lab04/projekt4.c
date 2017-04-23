#include <stdio.h>


enum kryterium
{
     MAX_F = 1,
     MAX_DELTA = 2,
     EUK_DELTA = 3
};
/******************************************************************************/

//------------------ Wartoœæ pierwszego równania x^2+y^2-8=0-------------------- 
double rownanie1(double x, double y)
{
       return x*x+y*y-8.0;
}

//------------------ Wartoœæ drugiego równania 4x+y^3=0-------------------- 
double rownanie2(double x, double y)
{
       return 4*x+y*y*y;;
}

/*============================= POCHODNE =====================================*/

//------------------ Wartoœæ pochodnej pierwszego równania wzglêdem x = 2*x -------------------- 
double rownanie1_dx(double x)
{
       return 2*x;
}

//------------------ Wartoœæ pochodnej drugiego równania wzglêdem x = 4 -------------------- 
double rownanie2_dx(double x)
{
       return 4;
}

//------------------ Wartoœæ pochodnej pierwszego równania wzglêdem y = 2*y -------------------- 
double rownanie1_dy(double y)
{
       return 2*y;
}

//------------------ Wartoœæ pochodnej drugiego równania wzglêdem y = 3y^2 -------------------- 
double rownanie2_dy(double y)
{
       return 3*y*y;
}

//------------------ Wartoœæ normy maksimum -------------------- 
double norma_max (double wektor[2])
{
       if(fabs(wektor[0]) > fabs(wektor[1]))
          return fabs(wektor[0]);
       
       return fabs(wektor[1]);
}

//------------------ Wartoœæ normy euklidesowej -------------------- 
double norma_euk(double wektor[2])
{
       return sqrt(wektor[0]*wektor[0]+wektor[1]*wektor[1]);
}

/*******************************************************************************
*******************************************************************************/

int main()
{
    int i=0;
    int kryt;
    double wyznacznik, eps;
    double x[2], x1[2], wart_fun[2], s[2];
    double macierz_poch[2][2], macierz_jacob[2][2];
    
    // ----- Podanie dok³adnoœci -------------
    printf("Podaj dokladnosc: ");
    scanf("%lf", &eps);
    
    // ----- Wybór kryterium zakoñczenia iteracji -------------
    printf("Podaj kryterium zakoñczenia iteracji: \n");
    printf(" 1 - ||f(x)|| max\n");
    printf(" 2 - ||x[n] - x[n-1]|| max\n");
    printf(" 3 - ||f(x)||_2 (Norma Euklidesowa)\n");
    scanf("%d", &kryt);
    
     // ----- Podanie wspó³rzêdnych punktu startowego -------------
    printf("Podaj wspolrzedne punktu startowego: \n");
    printf(" x = ");
    scanf("%lf", &x[0]);
    printf(" y = ");
    scanf("%lf", &x[1]);
    
    
    x1[0] = x[0] + 1;
    x1[1] = x[1] + 1;
    
    // ----- Obliczanie pierwiasta ukladu równañ metod¹ Newtona -------------
    
    for(i=0; i<200; i++)
    {
       x1[0]=x[0];
       x1[1]=x[1];
       
       //------ Macierz pochodnych ---------
       macierz_poch[0][0] = rownanie1_dx(x[0]);
       macierz_poch[0][1] = rownanie1_dy(x[1]);
       macierz_poch[1][0] = rownanie2_dx(x[0]);
       macierz_poch[1][1] = rownanie2_dy(x[1]);
       
       //------ Wektor wartoœci równañ ---------
       wart_fun[0] = rownanie1(x[0], x[1]);
       wart_fun[1] = rownanie2(x[0], x[1]);
       
       //------ Wartoœæ wyznacznika ---------
       wyznacznik = macierz_poch[0][0]*macierz_poch[1][1] - macierz_poch[1][0]*macierz_poch[0][1];
       
       if(fabs(wyznacznik) < 0.00000000000001)
       {
          printf("Brak rozwiazania, macierz osobliwa");
          system("pause");
          return -1;
       }
       
       //------ Poprawka ze wzgledu na x ---------
       s[0] = (macierz_poch[1][1]*(wart_fun[0])-(wart_fun[1])*macierz_poch[0][1])/wyznacznik;
       
       //------ Poprawka ze wzgledu na y ---------
       s[1] = (macierz_poch[0][0]*(wart_fun[1])-(wart_fun[0])*macierz_poch[1][0])/wyznacznik;
       
       
       printf("\nIteracja: %d\tx = %.10lf\ty = %.10lf", i, x1[0], x1[1]);
       
       //------ Nowe przyblizenie ---------
       x[0] = x1[0] -s[0];
       x[1] = x1[1] -s[1];
       
       
       if(kryt == MAX_F && norma_max(wart_fun) < eps)
               break;
       if(kryt == MAX_DELTA && norma_max(s) < eps)
               break;
       if(kryt == EUK_DELTA && norma_euk(s) < eps)
               break;
    }
    
    printf("\n\nWyniki: \n");
    printf(" x = %.10lf\n", x1[0]);
    printf(" y = %.10lf\n", x1[1]);
    printf(" Iteracje = %d\n", i);
    
    
    system("pause");
}

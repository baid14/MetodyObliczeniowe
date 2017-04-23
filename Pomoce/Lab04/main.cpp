#include <stdlib.h>
#include <stdio.h>
#include <math.h>
/*MAX liczba iteracji*/
#define MAX 1000

/*Dokladnosc obliczen */
#define DOKLADNOSC 0.000000000000001

void oblicz(double x, double y){

    int i=0;                   //licznik iteracji
    double rozwiazania[2];     //tablica rozwi¹zañ funkcji
    double poprawki[2];        //tablica "poprawek"
    double jakobian[2][2];     //macierz Jacobiego
    double wyznacznik;         //wyznacznik

    while(i<MAX){
        rozwiazania[0]= x*x+y*y-8;  //liczymy wartoœci obu funkcji w zadanych punktach
        rozwiazania[1]= y*y*y+4*x;  //liczymy wartoœci obu funkcji w zadanych punktach

        /*Sprawdzamy warunek zakoñczenia - je¿eli wartoœci funkcji bed¹ wystarczaj¹co blisko zera*/
        if(rozwiazania[0]>rozwiazania[1]){
            if(rozwiazania[0]<DOKLADNOSC){
                printf("\n\nKoniec po %d iteracjach\nz dokladnoscia %.15lf\n(rozw.[0]<DOKLADNOSC)\n",i,DOKLADNOSC);
                break;
            }

        }
        else
            if(rozwiazania[1]<DOKLADNOSC){
                printf("\n\nKoniec po %d iteracjach\nz dokladnoscia %.15lf\n(rozw.[1]<DOKLADNOSC)\n",i,DOKLADNOSC);
                break;
            }

        /*Tworzymy macierz Jacobiego - pochodne cz¹stkowe*/
        jakobian[0][0]=2*x;
        jakobian[0][1]=2*y;
        jakobian[1][0]=4;
        jakobian[1][1]=3*y*y;

        /*Liczymy i sprawdzamy wyznacznik - je¿eli jest równy zero to brak rozwi¹zañ*/
        wyznacznik = (jakobian[0][0]*jakobian[1][1]-jakobian[0][1]*jakobian[1][0]);
        if(wyznacznik==0){
            printf("\tMacierz jest osobliwa, brak rozwiazan!\n");
            break;
        }

        /*Liczymy kolejne "poprawki"*/
        poprawki[0] = (rozwiazania[0]*jakobian[1][1]-jakobian[0][1]*rozwiazania[1])/wyznacznik;
        poprawki[1] = (jakobian[0][0]*rozwiazania[1]-rozwiazania[0]*jakobian[1][0])/wyznacznik;

        /*Wnosimy poprawki do szukanych rozwi¹zañ*/
        x = x - poprawki[0];
        y = y - poprawki[1];

        /*Wypisujemy wyniki i liczbe iteracji*/
		printf("\nIteracja: %d \t",++i);
        printf("\n x: %.15lf\n y: %.15lf\n------------------------", x, y);
    }
}

int main(){
    printf("Metoda Newtona\n\n");
    oblicz((sqrt(8.0)), 0.001);
    return 0;
}

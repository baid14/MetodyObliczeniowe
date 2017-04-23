#include <cstdlib>
#include <iostream>
#include <cmath>
#include <stdio.h>
//#include "calerf.h"

using namespace std;

//deklaracja wskaznikow do plikow z danymi
    FILE *ogolne, *maxb, *maxbwykr, *wykr2, *wykr3;
//deklaracja stalych i zmiennych projektowych
    double h , dt, lambda, D=1;
    int x, t, a=1, b=11, tmax=2;

//tablica wynikowa dla laasonen
    double** wynik;
//tablica wynikowa dla crankaNicolson
    double** wynik2;

//tablica testowa
    double** testowe;

/** Algorytm Thomasa */
double *metodaThomasa(double **A, double *b,int m)
{
            double n,r,suma=0;
            // wektor wartosci wynikowych
            double* x = new double [m];
            // przeksztalcenie macierzy do postaci trojkatnej gornej
            for(int i=1;i<m;i++)
            {
                    n=A[i-1][i-1];
                    r=b[i-1];
                    A[i][i]=A[i][i]-A[i][i-1]*pow(n,-1)*A[i-1][i];
                    b[i]=b[i]-A[i][i-1]*pow(n,-1)*r;
                    A[i][i-1]=0;
            }
            // algorytm rozwiazywania ukladu row. liniowych z macierza trojkatna gorna
            x[m-1]=b[m-1]/A[m-1][m-1];
            for(int i=m-2;i>=0;i--)
            {
                    for(int j=i+1;j<m;j++)
                    {
                            suma=suma+A[i][j]*x[j];
                    }
                    x[i]=(b[i]-suma)/A[i][i];
                    suma=0;
            }
            return x;
}

/** Metoda uzupelniajaca tablice rozwiazaniami analitycznymi */
void rozwAnalit()
{
    double tt=dt;
    double dh;
    for(int i=1;i<t;i++){
            dh=1;
            for(int j=0;j<x;j++){
               testowe[i][j]=1.0-(((double)a/(double)dh)*erfc((double)(dh-a)/(2.0*sqrt(D*tt))));
               dh+=h;
            }
            tt+=dt;
    }
}

/** Metoda realizuj�ca algorytm lasonen */
void laasonen()
{
    // macierz trojdiagonalna wspolczynnikow przy niewiadomych (stala)
    double** A = new double *[x];
        for (int i=0;i<x;i++) A[i] = new double[x];
    // uzupelnienie macierzy zerami
    for(int i=0;i<x;i++) for(int j=0;j<x;j++) A[i][j]=0;

    // wektor zawierajacy znane wartosci rownania z poprzedniego poziomu czasowego
    double* b2 = new double [x];
    for(int i=0;i<x;i++) b2[i]=0;

    // wektor wynikowy ktory zostanie uzupeniony obliczonymi wartosciami rownania
    double* w = new double [x];
    for(int i=0;i<x;i++) w[i]=0;

    double dh;
    // petla uzupelniajaca macierz wspolczynnikow (zalezna od postaci rownania)
    for(int k=1;k<t;k++){

            A[0][0]=1; // d[0]
            A[x-1][x-1]=1;
            dh=h+1;
            for(int i=1;i<x-1;i++){
            // wspolczynniki:
            // dla roznicy progresywnej pierwszej pochodnej przestrzennej w rownaniu z problemu D
                    //A[i][i-1]=lambda;
                    //A[i][i]=-(1.0 + ( 2 + ( (2*h)/dh  ))*lambda);
                    //A[i][i+1]=(1 + ( (2*h)/dh  )) *lambda;

            // dla roznicy centralnej pierwszej pochodnej przestrzennej w rownaniu z problemu D
                    A[i][i-1]=lambda*(1.0 - (h/dh));
                    A[i][i]=-(1.0 + (2*lambda));
                    A[i][i+1]=lambda*(1.0 + (h/dh));
                    b2[i] = -wynik[k-1][i];
                    dh+=h;
            }
            b2[0] = wynik[k][0];
            b2[x-1] = wynik[k][x-1];

            // obliczenie wyniku dla aktualnego poziomu czasowego za pomoca metody Thomasa
            w=metodaThomasa(A,b2,x);

            for( int i = 1; i < x-1; i++ ) wynik[k][i] = w[i];
    }

    // zwolnienie pamieci
    for(int i=x-1;i>=0;i--) delete []A[i];
    delete []A;
    delete []b2;
    delete []w;
}

void crankaNicolson()
{
    // macierz trojdiagonalna wspolczynnikow przy niewiadomych (stala)
    // dla szukanego poziomu czasowego
    double** A = new double *[x];
        for (int i=0;i<x;i++) A[i] = new double[x];
    // uzupelnienie macierzy zerami
    for(int i=0;i<x;i++) for(int j=0;j<x;j++) A[i][j]=0;

    // macierz trojdiagonalna wspolczynnikow (stala)
    // dla poprzedniego (znanego) poziomu czasowego
    double** B = new double *[x];
        for (int i=0;i<x;i++) B[i] = new double[x];
    // uzupelnienie macierzy zerami
    for(int i=0;i<x;i++) for(int j=0;j<x;j++) B[i][j]=0;

    // wektor zawierajacy znane wartosci rownania z poprzedniego poziomu czasowego
    double* b2 = new double [x];
    for(int i=0;i<x;i++) b2[i]=0;

    // wektor bedacy suma wektora b2 z wartosciami warunkow brzegowych
    double* b3 = new double [x];
    for(int i=0;i<x;i++) b3[i]=0;

    // wektor wynikowy ktory zostanie uzupeniony obliczonymi wartosciami rownania
    double* w = new double [x];
    for(int i=0;i<x;i++) w[i]=0;

    double dh;
    // petla uzupelniajaca macierze wspolczynnikow (zalezna od postaci rownania)
    for(int k=1;k<t;k++){

            A[0][0]=1;
            A[x-1][x-1]=1;
            B[0][0]=0;
            B[x-1][x-1]=0;
            dh=h+1;
            for(int i=1;i<x-1;i++){
            // wspolczynniki:
            // dla roznicy progresywnej pierwszej pochodnej przestrzennej w rownaniu z problemu D
                    //A[i][i-1]=lambda/2.0;
                    //A[i][i]=-(1.0 + ( 1 + ( (2*h)/dh  ))*lambda);
                    //A[i][i+1]=(0.5 + ( (2*h)/dh  )) *lambda;

            // dla roznicy centralnej pierwszej pochodnej przestrzennej w rownaniu z problemu D
                    A[i][i-1]=lambda*(0.5 - (h/dh));
                    A[i][i]=-(1.0 + lambda);
                    A[i][i+1]=lambda*(0.5 + (h/dh));
                    b2[i] = wynik2[k-1][i];
                    dh+=h;
            }

            // uzupelnienie macierzy wspolczynnikow przy wartosciach znanych
            for(int i=1;i<x-1;i++){
                    B[i][i-1]=-lambda/2.0;
                    B[i][i]=-(1 - lambda);
                    B[i][i+1]=-lambda/2.0;
            }
            b2[0] = wynik2[k-1][0];
            b2[x-1] = wynik2[k-1][x-1];

            // mnozenie macierzy B z wektorem b2 oraz dodanie wartosci war. brzegowych
            double suma=0;
            for(int i=0;i<x;i++){
                    for(int j=0;j<x;j++) suma+=B[i][j]*b2[j];
                    b3[i]=suma;
                    suma=0;
            }
            b3[0]+=wynik2[k][0];
            b3[x-1]+=wynik2[k][x-1];

            // obliczenie wyniku dla aktualnego poziomu czasowego za pomoca metody Thomasa
            w=metodaThomasa(A,b3,x);

            for( int i = 1; i < x-1; i++ ) wynik2[k][i] = w[i];
    }
    // zwolnienie pamieci
    for(int i=x-1;i>=0;i--) delete []A[i];
    delete []A;
    for(int i=x-1;i>=0;i--) delete []B[i];
    delete []B;
    delete []b2;
    delete []b3;
    delete []w;
}

void wstep()
{
     cout<<"\n  Projekt zaliczeniowy laboratorium z przedmiotu 'Metody Obliczeniowe'\n"
         <<"  semestr letni 2011/2012\n"
         <<"  autor: Marcin Dzierwa\n"
         <<"  WFMiI Informatyka 2 rok\n"
         <<"  nr. albumu: 93978\n"
         <<"  9.06.2012\n\n"
         <<"Projekt D\n"
         <<"rownanie:  dU(x,t)/dt = D [ d2U(x,t)/dx2 + 2/x * dU(x,t)/dx ]"
         <<endl<<"-------------------------------------------------------------------------"<<endl;
}

int main(int argc, char *argv[])
{
    wstep();
    // krok przestrzenny
    h=0.1;

// pliki do ktorych zostana zapisane odpowiednie dane
    //maxb = fopen("bladmax_laasonen.txt","w");
    //maxb = fopen("bladmax_CN.txt","w");
    //maxbwykr = fopen("bladmax_laasonen.csv","w");
    //maxbwykr = fopen("bladmax_CN.csv","w");
    //wykr2 = fopen("wykres2laasonen.csv","w");
    //wykr2 = fopen("wykres2CN.csv","w");
    //wykr3 = fopen("wykres3laasonen.csv","w");
    //wykr3 = fopen("wykres3CN.csv","w");

//while(1){

    cout<<"\n1. test. Podaj h:"<<endl;
    //cin>>h;
    cout<<"\nlicze dla h = "<<h<<endl;
    cout<<"\nasd "<<sizeof(double)<<endl;
    dt=h*h;

// inicjalizacja zmiennych (krok dt obliczany tak aby lambda = 1 przy D = 1)
    x=(int)((b-a)/h)+1;
    t=(int)(tmax/dt)+1;
    lambda = D*(dt/(h*h));
    cout<<"\ndane x= "<<x<< " t " << t <<endl;

// tablica wynikowa laasonen
    wynik = new double *[t];
        for (int i=0;i<t;i++) wynik[i] = new double[x];
    for(int i=0;i<t;i++) for(int j=0;j<x;j++) wynik[i][j]=0;

// tablica wynikowa crankaNicolson
    wynik2 = new double *[t];
        for (int i=0;i<t;i++) wynik2[i] = new double[x];
    for(int i=0;i<t;i++) for(int j=0;j<x;j++) wynik2[i][j]=0;

// tablica testowa wartosci analitycznych
    testowe = new double *[t];
        for (int i=0;i<t;i++) testowe[i] = new double[x];
    for(int i=0;i<t;i++) for(int j=0;j<x;j++) testowe[i][j]=0;


// ---------------------------------------------------------------------------------
//uzupelnianie tablic o warunki brzegowe i poczatkowy
//warunek poczatkowy
    for(int i=0;i<x;i++) wynik[0][i]=1;
    for(int i=0;i<x;i++) wynik2[0][i]=1;

//1 warunek brzegowy
    for(int i=1;i<t;i++) wynik[i][0]=0;
    for(int i=1;i<t;i++) wynik2[i][0]=0;
    double tt=dt;

//2 warunek brzegowy
    for(int i=1;i<t;i++){
             wynik[i][x-1]=1.0-(((double)a/(double)b)*erfc(10.0/(2.0*sqrt(D*tt))));
             tt+=dt;
    }
    tt=dt;
    for(int i=1;i<t;i++){
             wynik2[i][x-1]=1.0-(((double)a/(double)b)*erfc(10.0/(2.0*sqrt(D*tt))));
             tt+=dt;
    }

// ---------------------------------------------------------------------------------
// uzupelnienie macierzy rozwiazan analitycznych
    rozwAnalit();

// ---------------------------------------------------------------------------------
// metody obliczenia danego rownania rozniczkowego czastkowego

    laasonen();
    crankaNicolson();

// ---------------------------------------------------------------------------------
// porownanie wynikow z wart. analitycznymi

    cout<<"\nlaasonen:\tcrank-nicolson:\tanalityczne:"<<endl;
    for(int i=0;i<x;i++) printf("%d %.20lf\t%.10lf\t%.10lf\n",i ,wynik[150][i],wynik2[2][i],testowe[150][i]);


// ---------------------------------------------------------------------------------
// wyliczenie bladow i zapis do plikow

//plik z wartosciami bladow dla wszystkich wartosci z tablicy

   /*ogolne = fopen("ogolne_CN.txt","w");

   fprintf(ogolne,"Wartosci i bledy dla metody Cranka-Nicolson przy kroku przestrzennym %lf i czasowym %lf\n\n",h,dt);
   double dtt=dt, dxx=1;
   for(int i=1;i<t;i++){
           dxx=1;
           fprintf(ogolne,"\nczas: %.4lf\n",dtt);
           for(int j=0;j<x;j++){
                   fprintf(ogolne,"x = %.4lf\t| obliczone: %.15lf\t| analityczne: %.15lf\t  | blad bezwzgl. = %.15lf\n",
                   dxx,wynik2[i][j],testowe[i][j],fabs(testowe[i][j]-wynik2[i][j]));
                   dxx+=h;
           }
           dtt+=dt;
   }
   fclose(ogolne);*/




//B��dy maxymalne dla tmax w funkcji kroku przestrzennego

    double maxblad=0;
   fprintf(maxb,"\nczas tmax: %d\n",tmax);
   //printf("\nczas tmax: %d\n",tmax);
   for(int j=0;j<x;j++){
           if(fabs(testowe[t-1][j]-wynik2[t-1][j])>maxblad) maxblad=fabs(testowe[t-1][j]-wynik2[t-1][j]);
           //dxx+=h;
   }
   fprintf(maxb,"h = %.4lf\tmax. blad bezwzgl. = %.20lf\n",h,maxblad);
    printf("h = %.4lf\tmax. blad bezwzgl. = %.20lf\n",h,maxblad);


  // double dxx=1, maxblad=0;
  /* for(int j=0;j<x;j++){
           if(fabs(testowe[t-1][j]-wynik2[t-1][j])>maxblad) maxblad=fabs(testowe[t-1][j]-wynik2[t-1][j]);
   }
   fprintf(maxbwykr,"%.3lf;%.20lf\n",h,maxblad);*/




//Wartosci analityczne i wyliczone dla kilku wartosci t

   /*double dxx=1;
   for(int j=0;j<x;j+=10){
       fprintf(wykr2,"%.4lf;%.15lf;%.15lf; ;%.15lf;%.15lf; ;%.15lf;%.15lf\n",
       dxx,testowe[1249][j],wynik2[1249][j],testowe[2499][j],wynik2[2499][j],testowe[3749][j],wynik2[3749][j]);
       dxx+=10*h;
   }*/




//B��dy maksymalne w funkcji czasu t

   /*double maxblad,dtt=dt;

   for(int i=1;i<t;i+=100){
       maxblad=0;
       for(int j=0;j<x;j++){
           if(fabs(testowe[i][j]-wynik[i][j])>maxblad) maxblad=fabs(testowe[i][j]-wynik[i][j]);
       }
       fprintf(wykr3,"%.4lf;%.15lf\n",dtt,maxblad);
       dtt+=100*dt;
   }*/


   /*double dxx=1, maxblad=0;
   for(int j=0;j<x;j++){
           if(fabs(testowe[t-1][j]-wynik2[t-1][j])>maxblad) maxblad=fabs(testowe[t-1][j]-wynik2[t-1][j]);
   }
   fprintf(maxbwykr,"%.3lf;%.20lf\n",h,maxblad);*/



// ---------------------------------------------------------------------------------
// czyszczenie pamieci

    for(int i=t-1;i>=0;i--) delete []wynik[i];
    delete []wynik;
    for(int i=t-1;i>=0;i--) delete []wynik2[i];
    delete []wynik2;
    for(int i=t-1;i>=0;i--) delete []testowe[i];
    delete []testowe;

    //if(h<=1&&h>0.11) h-=0.1;
    //else if(h<=0.11&&h>0.011) h-=0.01;
    //else if(h<=0.011) break;
//}

//zamykanie odpowiednich plikow

//fclose(maxb);
//fclose(maxbwykr);
//fclose(wykr2);
//fclose(wykr3);

    system("PAUSE");
    return EXIT_SUCCESS;
}

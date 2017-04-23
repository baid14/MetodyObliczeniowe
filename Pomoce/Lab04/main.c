	#include <math.h>
    #include <stdio.h>
    #include <stdlib.h>
    //dokladnosc
    const double eps = 0.000001;
    //ilosc obrotow
    const double loop = 50;
    double f1(double x, double y) {
            return x*x + y*y - 8;
    }
    double f2(double x, double y) {
            return 4*x + y*y*y;
    }
    void oblicz_pierwiastki(){
            int i = 1;
            double x[2] = {0, sqrt(8)};//przyblizenia poczatkowe
            double j[2][2];//macierz Jacobiego
            double tmp[2];//wektor pomocniczy/poprawek, wyliczany w kazdej iteracji
            double detJ, detJx, detJy, fx, fy;      
            printf("i\tx\t\ty\t\tf1(x,y)\t\tf2(x,y)\n");
            printf("*****************************************************************\n");
            while(1) {
                    //wartosc funkcji dla wyznaczonych x, y
                    fx = f1(x[0],x[1]);
                    fy = f2(x[0],x[1]);
                    //obliczenie macierzy Jacobiego z pochodnych
                    j[0][0] = 2*x[0];
                    j[0][1] = 2*x[1];                        
                    j[1][0] = 4;
                    j[1][1] = 3*x[1]*x[1];          
                    //wyznacznik macierzy Jacobiego
                    detJ = j[0][0] * j[1][1] - j[1][0] * j[0][1];
                    detJx = fx * j[1][1] - fy * j[0][1];
                    detJy = j[0][0] * fy - j[1][0] * fx;
                   
                    tmp[0] = detJx / detJ;
                    tmp[1] = detJy / detJ;
     
                    //kolejne przyblizenie
                    x[0] -= tmp[0];
                    x[1] -= tmp[1];
                    printf("%d\t%lf\t%lf\t%lf\t%lf\n",i,x[0],x[1],fx,fy);
                   
                    if(i >= loop){
                            printf("\nKONIEC OBLICZEN\nOsiagnieto maksymalna ilosc iteracji\n\n");
                            break;
                    }
                    if(fabs(fx) < eps && fabs(fy) < eps){
                            printf("\nKONIEC OBLICZEN\nOsiagnieto zadana dokladnosc dla wartosci funkcji\n\n");
                            break;
                    }
                    if(fabs(tmp[0]) < eps && fabs(tmp[1]) < eps){
                            printf("\nKONIEC OBLICZEN\nOsiagnieto zadana dokladnosc dla wektora poprawki\n\n");
                            break;
                    }
            i++;
            }
    }
    int main(void) {
            printf("Metoda Newtona\n");
            oblicz_pierwiastki();
            system("pause");
    }

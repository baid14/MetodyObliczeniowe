#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define H 20 // iteracje

double f1(double x, double y) {
	return x*x + y*y - 8;
}

double f2(double x, double y) {
	return 4*x + y*y*y;
}

//Macierz Jacobiego
int main(){
	double J[2][2], wek[2], w, fx1, fx2, p_x = 0, p_y = sqrt(8), residuum; //macierz Jacobiego, fx1, fx2 wartosci po policzeniu funkcji,
		   					   			 	   	  	  			 		   // w - wyznacznik macierzy, p_x, p_y - pocz¹tkoer przybli¿enia	
	int i=0;
	double eps=0.00001; // epsilon
	do {
		J[0][0] = 2*p_x; J[0][1] = 2*p_y;  J[1][0] = 4; J[1][1] = 3*p_y*p_y; //Jakobian	
		w = J[0][0] * J[1][1] - J[1][0] * J[0][1]; //obliczenie wyznacznika z macierzy Jakobiego
		fx1 = f1(p_x,p_y); //wartosc funkcji dla wyznaczonego x
		fx2 = f2(p_x,p_y); //wartosc funkcji dla wyznaczonego y
		//obliczanie wektora dla kolejnych iteracji
		wek[0]=(fx1 * J[1][1] - fx2 * J[0][1])/w;
		wek[1] = (J[0][0] * fx2 - J[1][0] * fx1)/w;
		//wyliczanie kolejnych przyblizeñ
		p_x -= wek[0];
		p_y -= wek[1];
		residuum = fmax(fabs(wek[0]), fabs(wek[1])); //wyliczenie residuum
		printf("Przyblizenie %d:\nx: %lf\ty: %lf\tresiduum: %lf\n****************************************************\n",i,p_x,p_y,residuum);
		//kryteria zakonczenia
		if(i>=H){
			printf("\nspelniony warunek liczby iteracji\n\n");
			break;
		}
		if(fabs(fx1)<eps && fabs(fx2)<eps){
			printf("\nspelniony warunek dokladnosci dla wartosci funkcji\n\n");
			break;
		}
		if(fabs(wek[0])<eps && fabs(wek[1])<eps){
			printf("\nspelniony warunek dokladnosci dla wektora iteracji\n\n");
			break;
		}
        i++;
	}while(1);

	system("PAUSE");
	return 0;
} 


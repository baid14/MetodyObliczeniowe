

#include <stdio.h>
#include <stdlib.h>

int main()
{
	//------- DANE - typ float
	float liczba_float = 1.0f;							// Liczba pojedynczej precyzji.
	float epsilon_float = liczba_float*0.5f;			// Wartosc poczatkowa epsilona float.
	float tmp_f = liczba_float + 1.0f;					// Zmienna pomocnicza dla float.
	float epsf_maszyn;									// Zmienna przechowujaca wartosc epsilona maszynowego float.

	//------- DANE - typ double
	double liczba_double = 1.0;							// Liczba podwojnej precyzji.
	double epsilon_double = liczba_double*0.5;			// Wartosc poczatkowa epsilona double.
	double tmp_d = liczba_double + 1.0;					// Zmienna pomocniacza dla double.
	double epsd_maszyn;									// Zmienna przechowujaca wartosc epsilona maszynowego double.

	int liczba_bitow = 0;								//Licznik bitow mantysy


//======================== FLOAT =====================================================================
	 
	/****************** Obliczanie liczby bitów mantysy i epsilona maszynowego - float *************/

	while(tmp_f > liczba_float)							// Petla bedzie sie wykonywaæ dopoki suma liczby i epsilona bedzie wieksza od tej liczby.
	{
		epsf_maszyn = epsilon_float;					// Zapamietanie wartosci epsilona maszynowego float.

		epsilon_float = epsilon_float*0.5f;				// Zmniejszenie wartosci epsilona float. Graniczna wartosc bledu.
		tmp_f = liczba_float + epsilon_float;			// Suma liczby i 'nowego' epsilona float.
		
		liczba_bitow++;									// Licznik bitow mantysy - dka float bedzie to 23.
	}
		

	/****************** Wydruk bitów mantysy i epsilona maszynowego - float *************/
	printf("------ FLOAT --------\n");
	printf(" Liczba bitow mantysy: %d\n Epsilon maszynowy: %g\n\n", liczba_bitow, epsf_maszyn);



//======================== DOUBLE =====================================================================

	/****************** Obliczanie liczby bitów mantysy i epsilona maszynowego - double *************/

	liczba_bitow = 0;									// Zerowanie licznika.

	while(tmp_d > liczba_double)						// Petla bedzie sie wykonywaæ dopoki suma liczby i epsilona bedzie wieksza od tej liczby.
	{
		epsd_maszyn = epsilon_double;					// Zapamietanie wartosci epsilona maszynowego double.

		epsilon_double = epsilon_double*0.5;			// Zmniejszenie wartosci epsilona double.
		tmp_d = liczba_double + epsilon_double;			// Suma liczby i 'nowego' epsilona double.

		liczba_bitow++;									// Licznik bitow mantysy - dla double bedzie to 52.
	}

	/****************** Wydruk bitów mantysy i epsilona maszynowego - double *************/

	printf("------ DOUBLE --------\n");
	printf("Liczba bitow mantysy: %d\nEpsilon maszynowy: %g\n\n", liczba_bitow, epsd_maszyn);


	system("pause");
	return 0;
}
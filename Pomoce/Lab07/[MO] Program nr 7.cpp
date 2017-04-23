#include <cstdlib>
#include <iostream>
#include <cmath>
#include <limits>
#include <iomanip>

#define ITERATIONS 64


namespace console
{
	void pause ()
	{
		std::cin.ignore (std::numeric_limits <std::streamsize>::max (), '\n');
		std::cout << "Nacisnij dowolny klawisz aby kontynuowac..." << std::endl;
		std::cin.get ();
	}
}


double b_VECTOR [2] = {49.0, 30.0};
double x_VECTOR [2] = {1.0, 1.0};

double compute_epsilon ()
{
    double epsilon = 1.0, temp = epsilon + 1.0;
    while (temp > 1.0)
    {
		epsilon /= 2.0;
		temp = epsilon + 1.0;
	}
	return epsilon;
}

double epsilon = 1e-6;//compute_epsilon ();


void fill_matrix (double * * matrix, int size)
{
	double MATRIX [4] = {5.0, 1.0, 4.0, 10.0};
	for (int i = 0; i < size; ++i)
		for (int j = 0; j < size; ++j)
			matrix [i] [j] = MATRIX [i * size + j];
	return;
}

void fill_vector (double * vector, double * source, int size)
{
	for (int i = 0; i < size; ++i)
		vector [i] = source [i];
	return;
}

void show_matrix (double * * matrix, int size)
{ 
	for (int i = 0; i < size; ++i)
	{
		for(int j = 0; j < size; ++j)
			std::cout << matrix [i] [j] << "\t";
		std::cout << std::endl;
	}
	std::cout << std::endl;
	return;
}

template <class type1>
void show_vector (type1 * vector, int size)
{
	for (int i = 0; i < size; ++i)
		std::cout << std::fixed << std::setprecision (14) << vector [i] << "\t";
	std::cout << std::endl << std::endl;
	return;
}

double compute_error (double * x, double * x_0, double * b, int size)
{
	double temp1, temp2 = temp1 = 0.0;
	for (int i = 0; i < size; ++i)
	{
		if (fabs (x_0 [i] - x [i]) > temp1)
			temp1 = fabs (x_0 [i] - x [i]);
		if (fabs (b [i]) > temp2)
			temp2 = fabs (b [i]);
	}
	return temp1 / temp2;
}

double * compute_residuum (double * * A, double * b, double * x, int size)
{
	// wektor residualnny
	// A * x_c - b = r
	
	double * residuum = new double [size];
	for (int i = 0; i < size; ++i)
	{
		residuum [i] = 0.0;
		for (int j = 0; j < size; ++j)
			residuum [i] += A [i] [j] * x [j];
		residuum [i] -= b [i];
	}
	return residuum;
}

void Jacobi_method (double * * A, double * b, double * x_0, double * x, int size)
{
	double sigma;
	for (int i = 0; i < size; ++i)
	{
		sigma = 0.0;
		for (int j = 0; j < i; ++j)
			sigma += (A [i] [j] * x [j]);
		
		for (int j = i + 1; j < size; ++j)
			sigma += (A [i] [j] * x [j]);
		
		x [i] = (b [i] - sigma) / A [i] [i];
	}
	return;
}

void Gauss_Seidel_method (double * * A, double * b, double * x_0, double * x, int size)
{
	double sigma;
	for (int i = 0; i < size; ++i)
	{
		sigma = 0.0;
		for (int j = 0; j < i; ++j)
			sigma += (A [i] [j] * x [j]);
		
		for (int j = i + 1; j < size; ++j)
			sigma += (A [i] [j] * x_0 [j]);
		
		x [i] = (b [i] - sigma) / A [i] [i];
	}
	return;
}

void SOR_method (double * * A, double * b, double * x_0, double * x, int size)
{
	double sigma, omega = 0.5;
	for (int i = 0; i < size; ++i)
	{
		sigma = 0.0;
		for (int j = 0; j < i; ++j)
			sigma += (A [i] [j] * x [j]);
		
		for (int j = i + 1; j < size; ++j)
			sigma += (A [i] [j] * x_0 [j]);
		
		sigma = (b [i] - sigma) / A [i] [i];
		x [i] = x_0 [i] + omega * (sigma - x_0 [i]);
	}
	return;
}

void iteration_method (void (* function) (double * * , double * , double * , double * , int ), int size)
{
	
	double * * A = new double * [size];
	for (int i = 0; i < size; ++i)
		A [i] = new double [size];
	fill_matrix (A, size);
	
	double * b = new double [size];
	fill_vector (b, b_VECTOR, size);
	
	double * x_0 = new double [size], * x = new double [size];
	
	for (int i = 0; i < size; ++i)
		x [i] = x_0 [i] = x_VECTOR [i];
	
	std::cout << "Wartosci poczatkowe wektora rozwiazan: " << std::endl;
	show_vector (x, size);
	
	int k = ITERATIONS;
	while (--k)                                                         // warunek wyjscia: liczba obiegow petli
	{
		function (A, b, x_0, x, size);
	
		if (compute_error (x, x_0, b, size) < epsilon)                  // warunek wyjscia: zakres bledu wiekszy niz tolerancja
			break;
		
		std::cout << ITERATIONS - k << ". krok" << std::endl
		          << "Wartosci wektora rozwiazan: " << std::endl;
		show_vector (x, size);
		std::cout << "Residuum ukladu: " << std::endl;
		show_vector (compute_residuum (A, b, x, size), size);
		
		for (int i = 0; i < size; ++i)
			x_0 [i] = x [i];
	}
	//print x
	
	delete [] x;
	delete [] x_0;
	delete [] b;
	for (int i = 0; i < size; ++i)
		delete [] A [i];
	delete [] A;
	return;
}

int main_menu ()
{
	int choice;
	std::cout << "\033[2J\033[1;1H";
	std::cout << "Rozwiazywanie ukladow rownan liniowych metodami iteracyjnymi" << std::endl << std::endl
	          << "Wybierz interesujaca Cie metode: " << std::endl
	          << "1. Metoda Jacobiego" << std::endl
	          << "2. Metoda Gaussa-Seidela" << std::endl
	          << "3. Metoda sukcesywnej nadrelaksacji" << std::endl
	          << "4. Wyjdz z programu" << std::endl << std::endl;			
	std::cin >> choice;	
	while (choice < 1 || choice > 4)
	{
		std::cout << "Podano zla wartosc, sprobuj ponownie" << std::endl;
		std::cin >> choice;
	}
	return choice;
}


int main (void)
{
	int size = 2;
	bool exit = false;
	while (!exit)
	{
		switch (main_menu ())
		{
			case 1:
				iteration_method (Jacobi_method, size);
				console::pause ();
				break;
			case 2:
				iteration_method (Gauss_Seidel_method, size);
				console::pause ();
				break;
			case 3:
				iteration_method (SOR_method, size);
				console::pause ();
				break;
			case 4:
				exit = true;
				break;
		}
	}
	system("PAUSE");
	return EXIT_SUCCESS;
}

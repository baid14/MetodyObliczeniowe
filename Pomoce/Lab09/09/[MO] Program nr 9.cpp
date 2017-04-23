#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <limits>
#include <iomanip>


namespace console
{
	void pause ()
	{
		std::cin.ignore (std::numeric_limits <std::streamsize>::max (), '\n');
		std::cout << "Nacisnij dowolny klawisz aby kontynuowac..." << std::endl;
		std::cin.get ();
	}
}


void gnuplot_script ()
{
	std::ofstream opened_file;
	opened_file.open (".data/plotter.gp");
	opened_file << "set title \"Dependency of the derivative error from a given step \\n in logarithmic scale\"" << std::endl
	            << "set encoding" << std::endl
	            << "set grid" << std::endl
	            << "set xrange [-15:0]" << std::endl
	            << "set yrange [-0.5:0]" << std::endl
	            << "set xlabel \"log_10 (h)\"" << std::endl
	            << "set ylabel \"log_10 (error)\"" << std::endl
	            << "set terminal png size 1024, 600" << std::endl
	            << "set output \"../.images/output.png\"" << std::endl
	            << "plot \"../.data/differential_method\" using 1:2 title \"Differential method\" with linespoints,\\" << std::endl
	            << "\"../.data/Numerov_method\" using 1:2 title \"Numerov method\" with linespoints" << std::endl;
	opened_file.close();
	return;
}

void print_vector (double * h_vector, double * vector, int size, std::string filename)
{
	std::ofstream opened_file;
	opened_file.open (filename.c_str ());
	opened_file << std::fixed << std::left << std::setprecision (20);
	for (int i = 0; i < size; ++i)
		opened_file << h_vector [i] << "\t" << vector [i] << std::endl;
	opened_file << std::endl;
	opened_file.close ();
	return;
}
void print_vector (double * h_vector, double * vector, int size)
{
	std::cout << std::fixed << std::setprecision (20);
	for (int i = 0; i < size; ++i)
		std::cout << h_vector [i] << "\t" << vector [i] << std::endl;
	std::cout << std::endl;
	std::cout.unsetf (std::ios_base::showbase);
	return;
}

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

int get_iterations ()
{
	int iterations = 1;
	double epsilon = compute_epsilon ();
	double h = 0.1;
	while (h > epsilon)
	{
		++iterations;
		h /= 10.0;
	}
	return iterations;
}

inline double analytic_method (double x)
{
	return (exp (x / 2.0) - exp (2.0 - x / 2.0)) / (1 - exp (2.0));
}

inline double *  get_analytic_solution (double a, double b, double h, int size)
{
	double * x = new double [size],
	       * analytic = new double [size];
	
	x [0] = a;
	x [size - 1] = b;
	
	analytic [0] = analytic_method (x [0]);
	for (int i = 1; i < size - 1; ++i)
	{
		x [i] = x [i - 1] + h;
		analytic [i] = analytic_method (x [i]);
	}
	analytic [size - 1] = analytic_method (x [size - 1]);
	
	delete [] x;
	return analytic;
}

//algorytm Thomasa
inline void tridiagonal_matrix_algorithm (double * l, double * d, double * u, double * x, double * b, int size)
{
	double * beta = new double [size];
	beta [0] = - u [0] / d [0];
	for (int i = 1; i < size; ++i)
		beta [i] = - u [i] / (l [i] * beta [i - 1] + d [i]);
	
	double * gamma = new double [size];
	gamma [0] = b [0] / d [0];
	for (int i = 1; i < size; ++i)
		gamma [i] = (b [i] - l [i] * gamma [i - 1]) / (l [i] * beta [i - 1] + d [i]);
	
	x [size - 1] = (b [size - 1] - l [size - 1] * gamma [size - 2]) / (l [size - 1] * beta [size - 2] + d [size - 1]) /*= gamma [size - 1]*/;
	for (int i = size - 2; i >= 0; --i)
		x [i] = beta [i] * x [i + 1] + gamma [i];
	
	delete [] gamma;
	delete [] beta;
	return;
}

inline double get_max_error (double * numerical, double a, double b, double h, int size)
{
	double max_error = 0.0,
	       * analytic = get_analytic_solution (a, b, h, size);
	for (int i = 0; i < size; ++i)
		max_error = std::max (max_error, fabs (analytic [i] - numerical [i]));
	delete [] analytic;
	delete [] numerical;
	return max_error;
}

double get_grade (double x_1, double x_2, double fx_1, double fx_2)
{
	return tan ((fx_1 - fx_2) / (x_1 - x_2));
}

// trzypunktowa dyskretyzacja konwencjonalna || metoda roznicowa
inline double * differential_method (double A, double B, double C, double h, double Y_a, double Y_b, int size)
{
	double * l = new double [size],
	       * d = new double [size],
	       * u = new double [size];
	
	u [size - 1] = l [0] = 0.0;
	d [size - 1] = d [0] = 1.0;
	
	double side_diagonal = A / (h * h), mid_diagonal = B - (2 * A / (h * h));
	for (int i = 1; i < size - 1; ++i)
		u [i - 1] = l [i] = side_diagonal;
	for (int i = 1; i < size - 1; ++i)
		d [i] = mid_diagonal;
	
	double * b = new double [size];
	b [0] = Y_a;
	b [size - 1] = Y_b;
	for (int i = 1; i < size - 1; ++i)
		b [i] = C;
	
	double * x = new double [size];
	
	tridiagonal_matrix_algorithm (l, d, u, x, b, size);
	
	delete [] b;
	delete [] u;
	delete [] d;
	delete [] l;
	
	return x;
}

// metoda Numerowa
inline double * Numerov_method (double A, double B, double C, double h, double Y_a, double Y_b, int size)
{
	double * l = new double [size];
	double * d = new double [size];
	double * u = new double [size];
	u [size - 1] = l [0] = 0.0;
	d [size - 1] = d [0] = 1.0;
	
	double mid_diagonal = (24.0 * A - 10.0 * B * h * h) / (12.0 * A - B * h * h);
	for (int i = 1; i < size - 1; ++i)
		u [i - 1] = l [i] = -1.0;
	for (int i = 1; i < size - 1; ++i)
		d [i] = mid_diagonal;
	
	double * b = new double [size];
	b [0] = Y_a;
	b [size - 1] = Y_b;
	for (int i = 1; i < size - 1; ++i)
		b [i] = C;
	double * x = new double [size];
	
	tridiagonal_matrix_algorithm (l, d, u, x, b, size);
	delete [] b;
	delete [] u;
	delete [] d;
	delete [] l;
	
	return x;
}

void numerical_method (double * (* function) (double , double , double , double , double , double , int ), std::string path)
{
	int h_steps = get_iterations (),                                // h_steps = 16
	    size;
	double a = 0.0, Y_a = 1.0,
	       b = 2.0, Y_b = 0.0,
	       h = 0.1,
	       A = 1.0, B = -0.25, C = 0.0;
	double * error_vector = new double [h_steps],
	       * h_vector = new double [h_steps];
	
	for (int i = 0; i < h_steps; ++i)
	{
		size = static_cast <int> ((fabs (b - a) / h) + 1.0);
		error_vector [i] = log10 (get_max_error (function (A, B, C, h, Y_a, Y_b, size), a, b, h, size));
		h_vector [i] = log10 (h);
		h /= 2.0;
	}
	
	std::cout << "Rzad dokladnosci wynosi: " << get_grade (h_vector [0], h_vector [1], error_vector [0], error_vector [1]) << std::endl;
	
	//print_vector (h_vector, error_vector, h_steps, path);
	print_vector (h_vector, error_vector, h_steps);
	
	delete [] error_vector;
	return;
}

int main_menu ()
{
	int choice;
	std::cout << "\033[2J\033[1;1H";
	std::cout << "Zastosowanie dwoch metod:" << std::endl
	          << "a) metody roznicowej," << std::endl
	          << "b) metody Numerowa," << std::endl
	          << "do rozwiazania rownania rozniczkowego zwyczajnego drugiego rzedu:" << std::endl
	          << "d^2 U (x)       1" << std::endl
	          << "---------   -   -   U (x)   =   0" << std::endl
	          << "  d x^2         4" << std::endl
	          << "0 <= x <= 2, U (0) = 1, U (2) = 0" << std::endl << std::endl;
	std::cout << "Wybierz metode:" << std::endl
	          << "1. Metoda roznicowa" << std::endl
	          << "2. Metoda Numerowa" << std::endl
	          << "3. Wyjdz z programu" << std::endl;
	std::cin >> choice;
	while (choice < 1 || choice > 3)
	{
		std::cout << "Podano zla wartosc, sprobuj ponownie" << std::endl;
		std::cin >> choice;
	}
	return choice;
}

int main (void)
{
	bool exit = false;
	do
	{
		switch (main_menu ())
		{
			case 1:
				numerical_method (differential_method, ".data/differential_method");
				console::pause ();
				break;
			case 2:
				numerical_method (Numerov_method, ".data/Numerov_method");
				console::pause ();
				break;
			case 3:
				exit = true;
				break;
		}
	}
	while (!exit);
	gnuplot_script ();
	return EXIT_SUCCESS;
}

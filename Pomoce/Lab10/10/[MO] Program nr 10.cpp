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
	opened_file << "set title \"Dependency of the maximum derivative error through given steps\\n in logarithmic scale\"" << std::endl
	            << "set encoding" << std::endl
	            << "set grid" << std::endl
	            << "set xrange [-15:0]" << std::endl
	            << "set yrange [-15:0]" << std::endl
	            << "set xlabel \"log_10 (h)\"" << std::endl
	            << "set ylabel \"log_10 (max_error)\"" << std::endl
	            << "set terminal png size 1024, 600" << std::endl
	            << "set output \"../.images/output.png\"" << std::endl
	            << "plot \"../.data/explicit_Euler_method\" using 1:2 title \"Explicit Euler method\" with linespoints,\\" << std::endl
	            << "\"../.data/implicit_Euler_method\" using 1:2 title \"Implicit Euler method\" with linespoints,\\" << std::endl
	            << "\"../.data/trapezoid_method\" using 1:2 title \"Trapezoid method\" with linespoints" << std::endl;
	opened_file.close();
	return;
}

void print_vectors (double * h_vector, double * vector, int size, std::string filename)
{
	std::ofstream opened_file;
	opened_file.open (filename.c_str ());
	opened_file << std::left << std::setprecision (20);
	for (int i = 0; i < size; ++i)
		opened_file << h_vector [i] << "\t" << vector [i] << std::endl;
	opened_file << std::endl;
	opened_file.close ();
	return;
}

void print_vectors (double * h_vector, double * vector, int size)
{
	std::cout << std::left << std::setprecision (25);
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
		h /= 2.0;
	}
	return iterations;
}

inline double analytic_method (double t)
{
	return (1.0 + t) * exp (-t);
}

inline double f (double y_n, double t_n)
{
	return exp (-t_n) - y_n;
}

// bezposrednia metoda Eulera || explicit Euler method || forward Euler method
// y_{n + 1} = y_n + h * f (y_n, t_n)
inline double explicit_Euler_method (double y_n, double t_n, double dt)
{
	return y_n + dt * f (y_n, t_n);
}

// posrednia metoda Eulera || implicit Euler method || backward Euler method
// y_{n + 1} = y_n + h * f (y_{n + 1}, t_{n + 1})
inline double implicit_Euler_method (double y_n, double t_n, double dt)
{
	return (y_n + dt * exp (-(t_n + dt))) / (1.0 + dt);
}

// metoda trapezow || improved/modified Euler method || explicit trapezoidal rule || Heun's method
// y_{n + 1} = y_n + {h over 2} * (f (y_n, t_n) + f (y_{n + 1}, t_{n + 1}))
inline double trapezoid_method (double y_n, double t_n, double dt)
{
	return (y_n + (dt / 2.0) * (f (y_n, t_n) + exp (-(t_n + dt)))) / (1.0 - (dt / 2.0));
}

double get_grade (double x_1, double x_2, double fx_1, double fx_2)
{
	return tan ((fx_1 - fx_2) / (x_1 - x_2));
}

void numerical_method (double (* function) (double, double , double ), std::string path)
{
	int h_steps = get_iterations ();
	int t_steps = 100;
	double y_n = 0.0, dy,
	       t_n, dt = 0.1;
	
	double max_t;
	double * error_vector = new double [h_steps];
	double * h_vector = new double [h_steps];
	double * each_t = new double [t_steps];
	
	for (int i = 0; i < h_steps; ++i)
	{
		t_n = 0.0;
		max_t = 0.0;
		for (int j = 0; j < t_steps; ++j)
		{
			max_t = std::max (max_t, fabs (analytic_method (t_n) - function (y_n, t_n, dt)));
			t_n += dt;
		}
		error_vector [i] = log10 (max_t);
		h_vector [i] = log10 (dt);
		dy = f (y_n, t_n);
		y_n += dy;
		dt /= 2.0;
	}
	
	std::cout << "Rzad dokladnosci wynosi: " << get_grade (h_vector [0], h_vector [1], error_vector [0], error_vector [1]) << std::endl;
	
	//print_vectors (h_vector, error_vector, h_steps, path);
	print_vectors (h_vector, error_vector, h_steps);
	
	delete [] each_t;
	delete [] h_vector;
	delete [] error_vector;
	
	return;
}

void demonstrate_instability (double (* function) (double, double, double ), std::string path)
{
	int t_steps = 100;
	double y_n = 0.0, dy,
	       t_n, dt = 1.9;
	
	double * error_vector = new double [t_steps];
	double * t_vector = new double [t_steps];
	
	t_n = 0.0;
	for (int i = 0; i < t_steps; ++i)
	{
		error_vector [i] = fabs (analytic_method (t_n) - function (y_n, t_n, dt));
		t_vector [i] = t_n;
		dy = f (y_n, t_n);
		y_n += dy;
		t_n += dt;
	}
	
	print_vectors (t_vector, error_vector, t_steps, path + "1");
	
	
	
	y_n = 0.0;
	dt = 2.1;
	t_n = 0.0;
	for (int i = 0; i < t_steps; ++i)
	{
		error_vector [i] = fabs (analytic_method (t_n) - function (y_n, t_n, dt));
		t_vector [i] = t_n;
		dy = f (y_n, t_n);
		y_n += dy;
		t_n += dt;
	}
	
	print_vectors (t_vector, error_vector, t_steps, path + "2");
	
	delete [] t_vector;
	delete [] error_vector;
	
	return;
}

int main_menu ()
{
	int choice;
	std::cout << "\033[2J\033[1;1H";
	std::cout << "Zastosowanie trzech metod:" << std::endl
	          << "a) bezposredniej Eulera," << std::endl
	          << "b) posredniej Eulera," << std::endl
	          << "c) metody trapezow," << std::endl
	          << "do rozwiazania rownania rozniczkowego zwyczajnego pierwszego rzedu:" << std::endl
	          << "dy (t)" << std::endl
	          << "------   +   y (t)   -   exp (-t)   =   0" << std::endl
	          << "  dt  " << std::endl
	          << "t >= 0, y (0) = 1" << std::endl << std::endl;
	std::cout << "Wybierz metode:" << std::endl
	          << "1. Bezposrednia metoda Eulera" << std::endl
	          << "2. Posrednia metoda Eulera" << std::endl
	          << "3. Metoda trapezow" << std::endl
	          << "4. Demonstracja warunkowej stabilnosci bezposredniej metody Eulera" << std::endl
	          << "5. Wyjdz z programu" << std::endl;
	std::cin >> choice;
	while (choice < 1 || choice > 5)
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
				numerical_method (explicit_Euler_method, ".data/explicit_Euler_method");
				break;
			case 2:
				numerical_method (implicit_Euler_method, ".data/implicit_Euler_method");
				break;
			case 3:
				numerical_method (trapezoid_method, ".data/trapezoid_method");
				break;
			case 4:
				demonstrate_instability (explicit_Euler_method, ".data/instability");
				break;
			case 5:
				exit = true;
				break;
		}
		console::pause ();
	}
	while (!exit);
	gnuplot_script ();
	system("PAUSE");
	return EXIT_SUCCESS;
}

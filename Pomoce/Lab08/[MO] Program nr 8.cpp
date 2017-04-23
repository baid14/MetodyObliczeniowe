#define _USE_MATH_DEFINES
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


template <class type1>
inline type1 function (type1 x)
{
	return cos (x);
}

template <class type1>
inline type1 derivative_analytic (type1 x)
{
	return -sin (x);
}

template <class type1>
inline type1 difference_forward_bipoint (type1 x_1, type1 x_2)
{
	return (function (x_2) - function (x_1)) / (x_2 - x_1);
}

template <class type1>
inline type1 difference_forward_tripoint (type1 x_0, type1 x_1, type1 x_2)
{
	return (-function (x_2) + 4.0 * function (x_1) - 3.0 * function (x_0)) / (x_2 - x_0);
}

template <class type1>
inline type1 difference_backward_bipoint (type1 x_0, type1 x_1)
{
	return (function (x_1) - function (x_0)) / (x_1 - x_0);
}

template <class type1>
inline type1 difference_backward_tripoint (type1 x_0, type1 x_1, type1 x_2)
{
	return (3.0 * function (x_2) - 4.0 * function (x_1) + function (x_0)) / (x_2 - x_0);
}

template <class type1>
inline type1 difference_central_bipoint (type1 x_0, type1 x_2)
{
	return (function (x_2) - function (x_0)) / (x_2 - x_0);
}


template <class type1>
inline type1 logarithm (type1 x)
{
	return (fabs (x) > 0.0) ? log10l (x) : 0.0l;
}

template <class type1>
type1 compute_epsilon ()
{
    type1 epsilon = static_cast <type1> (1.0), temp = epsilon + static_cast <type1> (1.0);
    while (temp > static_cast <type1> (1.0))
    {
		epsilon /= static_cast <type1> (2.0);
		temp = epsilon + static_cast <type1> (1.0);
	}
	return epsilon;
}

template <class type1>
int get_iterations ()
{
	int iterations = 1;
	type1 epsilon = compute_epsilon <type1> ();
	type1 h = static_cast <type1> (0.1);
	while (h > epsilon)
	{
		++iterations;
		h /= static_cast <type1> (10.0);
	}
	return iterations;
}


template <class type1>
void get_relative_error_table (std::string filename)
{
	type1 a = 0.0;
	type1 b = M_PI_4; //M_PI_4; LUB                             atan (static_cast <type1> (1.0)); LUB atan2l (static_cast <type1> (0.0), static_cast <type1> (-1.0));
	type1 c = M_PI_2; //M_PI_2; LUB static_cast <type1> (2.0) * atan (static_cast <type1> (1.0)); LUB atan2l (static_cast <type1> (0.0), static_cast <type1> (-1.0));
	type1 h = 0.1;
	int size = get_iterations <type1> ();
	type1 * * relative_error_table = new type1 * [10];
	for (int i = 0; i < 10; ++i)
		relative_error_table [i] = new type1 [size];
	for (int i = 0; i < size; ++i)
	{
		relative_error_table [0] [i] = log10l (h);
		                                                                                                               // INTERVAL STARTPOINT
		relative_error_table [1] [i] = log10l (fabsl (derivative_analytic (a) - difference_forward_bipoint   (a, a + h)));                // forward difference,  2-point
		relative_error_table [2] [i] = log10l (fabsl (derivative_analytic (a) - difference_forward_tripoint  (a, a + h, a + static_cast <type1> (2.0) * h)));   // forward difference,  3-point
		
		                                                                                                               // INTERVAL CENTERPOINT
		relative_error_table [3] [i] = log10l (fabsl (derivative_analytic (b) - difference_backward_bipoint  (b, b + h)));                // backward difference, 2-point
		relative_error_table [4] [i] = log10l (fabsl (derivative_analytic (b) - difference_central_bipoint   (b - h, b + h)));            // central difference,  2-point
		relative_error_table [5] [i] = log10l (fabsl (derivative_analytic (b) - difference_forward_bipoint   (b, b + h)));                // forward difference,  2-point
        relative_error_table [6] [i] = log10l (fabsl (derivative_analytic (b) - difference_backward_tripoint (b - static_cast <type1> (2.0) * h, b - h, b)));   // backward difference, 3-point
		relative_error_table [7] [i] = log10l (fabsl (derivative_analytic (b) - difference_forward_tripoint  (b - static_cast <type1> (2.0) * h, b - h, b)));   // forward difference,  3-point
		
		                                                                                                               // INTERVAL ENDPOINT
		relative_error_table [8] [i] = log10l (fabsl (derivative_analytic (c) - difference_backward_bipoint  (c - h, c)));                  // backward difference, 2-point
		relative_error_table [9] [i] = log10l (fabsl (derivative_analytic (c) - difference_backward_tripoint (c - static_cast <type1> (2.0) * h, c - h, c)));    // backward difference, 3-point
		h /= static_cast <type1> (10.0);
	}
	print_relative_error_table (relative_error_table, size, filename);
	
	std::cout << "Kat nachylenia:" << std::endl;
	type1 angle;
	
	for (int i = 1; i < 9; ++i)
	{
		angle = fabsl (tan (relative_error_table [i] [0] - relative_error_table [i] [1]) / (relative_error_table [0] [0] - relative_error_table [0] [1]));
		std::cout << angle << std::endl;
	}
	
	for (int i = 0; i < 10; ++i)
		delete [] relative_error_table [i];
	delete [] relative_error_table;
	return;
}

template <class type1>
void print_relative_error_table (type1 * * relative_error_table, int size, std::string filename)
{
	std::string names [] = {"h", "forward difference", "forward difference", "backward difference", "central difference", "forward difference", "backward difference", "forward difference",  "backward difference", "backward difference"};
	std::string placement [] = {"h", "interval startpoint", "interval startpoint", "interval centerpoint", "interval centerpoint", "interval centerpoint", "interval centerpoint", "interval centerpoint", "interval endpoint", "interval endpoint"};
	std::string approximation [] = {"h", "bipoint approximation", "tripoint approximation", "bipoint approximation", "bipoint approximation", "bipoint approximation", "tripoint approximation", "tripoint approximation", "bipoint approximation", "tripoint approximation"};
	
	for (int i = 0; i < 10; ++i)
		std::cout << placement [i] << ": " << names [i] << "(" << approximation [i] << ")" << "\t";
	std::cout << std::endl;
	for (int j = 0; j < size; ++j)
	{
		for (int i = 0; i < 10; ++i)
			std::cout << std::fixed << std::setw (10) << std::setprecision (8) << relative_error_table [i] [j] << "\t";
		std::cout << std::endl;
	}
	
	std::ofstream opened_file;
	opened_file.open (("data/" + filename).c_str ());
	for (int i = 0; i < 10; ++i)
		opened_file << placement [i] << ": " << names [i] << "(" << approximation [i] << ")" << "\t";
	opened_file << std::endl;
	for (int j = 0; j < size; ++j)
	{
		for (int i = 0; i < 10; ++i)
			opened_file << std::fixed << std::setw (15) << std::setprecision (15) << relative_error_table [i] [j] << "\t";
		opened_file << std::endl;
	}
	opened_file.close ();
	
	
	opened_file.open (("gnuplot_files/" + filename + "_interval_startpoint").c_str ());
	opened_file << "set title \"Dependency of the derivative error from a given step \\n on the interval startpoint" << std::endl
	            << "set encoding" << std::endl
	            << "set xrange [-15.0:0.0]" << std::endl
	            << "set yrange [-12.0:1.0]" << std::endl
	            << "set xlabel \"log_10 (h)\"" << std::endl
	            << "set ylabel \"log_10 (error)\"" << std::endl
	            << "set grid" << std::endl
	            << "set output \"../images/" << filename << "_interval_startpoint.png\"" << std::endl
	            << "set terminal png" << std::endl
	            << "plot \"../data/" << filename << "\" using 1:2 title \"Forward 2P\" with linespoints,\\" << std::endl
	            << "\"../data/" << filename << "\" using 1:3 title \"Forward 3P\" with linespoints" << std::endl;
	opened_file.close ();
	
	opened_file.open (("gnuplot_files/" + filename + "_interval_centerpoint").c_str ());
	opened_file << "set title \"Dependency of the derivative error from a given step \\n on the interval startpoint" << std::endl
	            << "set encoding" << std::endl
	            << "set xrange [-15.0:0.0]" << std::endl
	            << "set yrange [-12.0:1.0]" << std::endl
	            << "set xlabel \"log_10 (h)\"" << std::endl
	            << "set ylabel \"log_10 (error)\"" << std::endl
	            << "set grid" << std::endl
	            << "set output \"../images/" << filename << "_interval_centerpoint.png\"" << std::endl
	            << "set terminal png" << std::endl
	            << "plot \"../data/" << filename << "\" using 1:4 title \"Backward 2P\" with linespoints,\\" << std::endl
	            << "\"../data/" << filename << "\" using 1:5 title \"Central 2P\" with linespoints,\\" << std::endl
	            << "\"../data/" << filename << "\" using 1:6 title \"Forward 2P\" with linespoints,\\" << std::endl
	            << "\"../data/" << filename << "\" using 1:7 title \"Backward 3P\" with linespoints,\\" << std::endl
	            << "\"../data/" << filename << "\" using 1:8 title \"Forward 3P\" with linespoints" << std::endl;
	opened_file.close ();
	
	opened_file.open (("gnuplot_files/" + filename + "_interval_endpoint").c_str ());
	opened_file << "set title \"Dependency of the derivative error from a given step \\n on the interval startpoint" << std::endl
	            << "set encoding" << std::endl
	            << "set xrange [-15.0:0.0]" << std::endl
	            << "set yrange [-12.0:1.0]" << std::endl
	            << "set xlabel \"log_10 (h)\"" << std::endl
	            << "set ylabel \"log_10 (error)\"" << std::endl
	            << "set grid" << std::endl
	            << "set output \"../images/" << filename << "_interval_endpoint.png\"" << std::endl
	            << "set terminal png" << std::endl
	            << "plot \"../data/" << filename << "\" using 1:9 title \"Backward 2P\" with linespoints,\\" << std::endl
	            << "\"../data/" << filename << "\" using 1:10 title \"Backward 3P\" with linespoints" << std::endl;
	opened_file.close ();
	return;
}

int first_menu ()
{
	int choice;
	std::cout << "\033[2J\033[1;1H";
	std::cout << "Wyznaczanie wartosci funkcji przy uzyciu wzorow roznicowych" << std::endl << std::endl
	          << "Wybierz typ zmiennych:" << std::endl
	          << "1. float" << std::endl
	          << "2. double" << std::endl
	          << "3. long double" << std::endl
	          << "4. Wyjdz z programu" << std::endl;
	std::cin >> choice;
	while (choice < 1 || choice > 4)
	{
		std::cout << "Wybrano niewlasciwy numer, sprobuj ponownie" << std::endl;
		std::cin >> choice;
	}
	return choice;
}

int main (void)
{
	bool exit = false;
	while (!exit)
	{
		switch (first_menu ())
		{
			case 1:
				get_relative_error_table <float> ("float");
				break;
			case 2:
				get_relative_error_table <double> ("double");
				break;
			case 3:
				get_relative_error_table <long double> ("long double");
				break;
			case 4:
				exit = true;
				break;
		}
		console::pause ();
	}
	return EXIT_SUCCESS;
}

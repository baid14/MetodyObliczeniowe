#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

template<class T> T horner(T x, int* a, const int n) {
		T y = static_cast<T>(a[0]);
		for(int i = 1; i < n; ++i) {
			y = y * x + a[i];
		}
		return y;
}
template<class T>T blad(T real, T y) {
		return static_cast<T>(fabs(real - y));
}
int main() {
		int a[] = {1,2,-3,4,-5};
    	int b[] = {8118, -11482, 1, 5741, -2030};
		const int n=5;
		double x1=2.0;
    	double r1=23.0;
		float x2=2.0;
    	float r2=23.0;
		long double x3=2.0;
    	long double r3=23.0;


		cout.precision(20);
		cout << "Wielomian:\ny = x^4 + 2x^3 - 3x^2 + 4x - 5\n\n";
		cout << "double:\nx = 2\ty = " << horner(x1, a, n) << "\tblad = "<< blad(r1, horner(x1, a, n)) << "\n\n"
			 << "float:\nx = 2\ty = " << horner(x2, a, n) << "\tblad = " << blad(r2, horner(x2, a, n)) << "\n\n"
			 << "long double:\nx = 2\ty = "<< horner(x3, a, n) << "\tblad = " << blad(r3, horner(x3, a, n))<<"\n"
			 << "\n\n";

		x1 = 0.707107;
		x2 = 0.707107;
		x3 = 0.707107;
		r1 = -1.9152732527e-07;
		r2 = -1.9152732527e-07;
		r3 = -1.9152732527e-07;

		cout.precision(20);
		cout << "Wielomian:\ny = 8118x^4 - 11482x^3 + x^2 + 5741x - 2030\n\n";
		cout << "double:\nx = 0.707107\ty = " <<horner(x1, b, n) << "\tblad = "<< blad(r1, horner(x1, b, n)) << "\n\n"
			 << "float:\nx = 0.707107\ty = " << horner(x2, b, n) << "\tblad = " << blad(r2, horner(x2, b, n)) << "\n\n"
			 << "long double:\nx = 0.707107\ty =  "<< horner(x3, b, n) << "\tblad = " << blad(r3, horner(x3, b, n)) << "\n\n";

		//system("PAUSE");
		return 0;
}


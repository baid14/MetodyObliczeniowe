#include <iostream>
#include <cmath>
#include<iomanip>

#define NMAX 50
#define TOLX 0.000001
#define TOLF 0.000001

using namespace std;

double rownanie1(double x, double y, double z)
{
	return x*x + y*y - 1.0;
}

double rownanie2(double x, double y, double z)
{
	return 2.0*x*x - y + 0.5;
}

double rownanie3(double x, double y, double z)
{
	return z*z*z*z*z - tan(x/4.0) - y*y*y;
}

double rownanie1pox(double x, double y, double z)
{
	return 2.0*x;
}

double rownanie1poy(double x, double y, double z)
{
	return 2.0*y;
}

double rownanie1poz(double x, double y, double z)
{
	return 0.0;
}

double rownanie2pox(double x, double y, double z)
{
	return 4.0*x;
}

double rownanie2poy(double x, double y, double z)
{
	return -1.0;
}

double rownanie2poz(double x, double y, double z)
{
	return 0.0;
}

double rownanie3pox(double x, double y, double z)
{
	return -0.25 * (1.0 + tan(x/4.0)*tan(x/4.0));
}

double rownanie3poy(double x, double y, double z)
{
	return -3.0*y*y;
}

double rownanie3poz(double x, double y, double z)
{
	return 5.0*z*z*z*z;
}

double nmax(double wektor[3])
{
	if(fabs(wektor[0]) > fabs(wektor[1]) && fabs(wektor[0]) > fabs(wektor[2]))
		return fabs(wektor[0]);

	if(fabs(wektor[1]) > fabs(wektor[0]) && fabs(wektor[1]) > fabs(wektor[2]))
		return fabs(wektor[1]);

	return fabs(wektor[2]);
}

int main()
{
 double xn[3] = {1.0, 1.0, 1.0},
        xn1[3],
        jacobi[3][3],
        fodx[3],
        delta[3],
        wyz_glowny,wyz_x,wyz_y,wyz_z;

     int i = 0;
    cout<<setw(10)<<"n"<<setw(15)<<"x1"<<setw(15)<<"x2"<<setw(15)<<"x3"<<setw(15)<<"EST(xn-x*)"<<setw(15)<<"|f(xn)|"<<endl<<endl;
    do
	{
        jacobi[0][0] = rownanie1pox(xn[0], xn[1], xn[2]);
		jacobi[0][1] = rownanie1poy(xn[0], xn[1], xn[2]);
		jacobi[0][2] = rownanie1poz(xn[0], xn[1], xn[2]);

		jacobi[1][0] = rownanie2pox(xn[0], xn[1], xn[2]);
		jacobi[1][1] = rownanie2poy(xn[0], xn[1], xn[2]);
		jacobi[1][2] = rownanie2poz(xn[0], xn[1], xn[2]);

		jacobi[2][0] = rownanie3pox(xn[0], xn[1], xn[2]);
		jacobi[2][1] = rownanie3poy(xn[0], xn[1], xn[2]);
		jacobi[2][2] = rownanie3poz(xn[0], xn[1], xn[2]);

        fodx[0] = -rownanie1(xn[0],xn[1],xn[2]);
        fodx[1] = -rownanie2(xn[0],xn[1],xn[2]);
        fodx[2] = -rownanie3(xn[0],xn[1],xn[2]);

        wyz_glowny =
             jacobi[0][0] * jacobi[1][1] * jacobi[2][2] + jacobi[0][1] * jacobi[1][2] * jacobi[2][0] + jacobi[0][2] * jacobi[1][0] * jacobi[2][1]
         -  (jacobi[2][0] * jacobi[1][1] * jacobi[0][2] + jacobi[2][1] * jacobi[1][2] * jacobi[0][0] + jacobi[2][2] * jacobi[1][0] * jacobi[0][1]);
cout <<"glowny  " <<wyz_glowny <<endl;
		wyz_x =
           fodx[0] * jacobi[1][1] * jacobi[2][2] + jacobi[0][1] * jacobi[1][2] * fodx[2] + jacobi[0][2] * fodx[1] * jacobi[2][1]
		- (fodx[2] * jacobi[1][1] * jacobi[0][2] + jacobi[2][1] * jacobi[1][2] * fodx[0] + jacobi[2][2] * fodx[1] * jacobi[0][1]);
cout <<"x   " <<wyz_x <<endl;
for ( int j = 0; j < 3; ++j )
		{
			for ( int k = 0; k < 3; k++ )
			{
				cout << jacobi[j][k] << " ";
			}cout << endl;
		}
		wyz_y =
           jacobi[0][0] * fodx[1] * jacobi[2][2] + fodx[0] * jacobi[1][2] * jacobi[2][0] + jacobi[0][2] * jacobi[1][0] * fodx[2]
		- (jacobi[2][0] * fodx[1] * jacobi[0][2] + fodx[2] * jacobi[1][2] * jacobi[0][0] + jacobi[2][2] * jacobi[1][0] * fodx[0]);
cout <<"y   " <<wyz_y <<endl;
		wyz_z =
           jacobi[0][0] * jacobi[1][1] * fodx[2] + jacobi[0][1] * fodx[1] * jacobi[2][0] + fodx[0] * jacobi[1][0] * jacobi[2][1]
		- (jacobi[2][0] * jacobi[1][1] * fodx[0] + jacobi[2][1] * fodx[1] * jacobi[0][0] + fodx[2] * jacobi[1][0] * jacobi[0][1]);
cout <<"z   " <<wyz_y <<endl;
		delta[0]= wyz_x/wyz_glowny;
		delta[1]= wyz_y/wyz_glowny;
		delta[2]= wyz_z/wyz_glowny;

		xn1[0] = xn[0] + delta[0];
		xn1[1] = xn[1] + delta[1];
		xn1[2] = xn[2] + delta[2];

		xn[0] = xn1[0];
		xn[1] = xn1[1];
		xn[2] = xn1[2];
        i++;

        //cout<<setw(10)<<i<<setw(15)<<xn[0]<<setw(15)<<xn[1]<<setw(15)<<xn[2]<<setw(15)<<nmax(delta)<<setw(15)<<nmax(fodx)<<endl;

	}while (((nmax(delta) > TOLX )|| (nmax(fodx) > TOLF )) && i<NMAX);

	return 0;
}

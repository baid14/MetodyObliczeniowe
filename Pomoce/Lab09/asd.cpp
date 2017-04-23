#include <iostream>
#include <cmath>
#include <cstdio>
#define M_PI 3.14159265358979323846

using namespace std;

void eliminacja(double *l, double *u, double *d, double *tmp_d, int n)
{
    for(int i=1;i<n;i++)
    {
        tmp_d[i]=d[i]-(l[i-1]*u[i-1]/tmp_d[i-1]);
    }
}
void thomas(double *l, double *u, double *d,double *b, double *x, int n)
{
    double *tmp_d = new double[n];
    double *tmp_b = new double[n];
    tmp_d[0]=d[0];
    tmp_b[0]=b[0];
    eliminacja(l,u,d,tmp_d,n);
    for(int i=1;i<n;i++)
    {
        tmp_b[i]=b[i]-(l[i-1]*tmp_b[i-1]/tmp_d[i-1]);
    }
    x[n-1]=tmp_b[n-1]/tmp_d[n-1];
    for(int i=n-2;i>=0;i--)
    {
        x[i]=(tmp_b[i]-u[i]*x[i+1])/tmp_d[i];
    }
    delete[] tmp_d;
    delete[] tmp_b;
}

void dyskretyzacja_numerowa(double h, double y0,double yn, int N, double *l, double *u, double *d,double *b, double *x)
{
    d[0] = 1.0;
    u[0] = 0.0;
    b[0] = y0;
    double h_kw=h*h;
    for(int i=1;i<N-1;i++)
    {
        l[i-1]=1.0-h_kw/12.0;
        d[i]=-(2+10.0*h_kw/12.0);
        u[i]=1.0-h_kw/12.0;
        b[i]=(sin((i-1)*h)+sin(i*h)*h_kw*10.0+sin((i+i)*h))*h_kw/6.0;
    }
    d[N-1]=1.0;
    l[N-2]=0.0;
    b[N-1]=yn;
    thomas(l,u,d,b,x,N);
}


void dyskretyzacja_konwencjonalna(double h, double y0,double yn, int N, double *l, double *u, double *d,double *b, double *x)
{
    //l[0] = x0;
    //double h=(xn-x0)/N;
    double alfa = 0.0;
    double beta = 1.0;
    double gamma = -y0;
    cout << h<<endl;

    double fi = 0.0;
    double ni = 1.0;
    double theta = -yn;

    double p = 1.0;
    double q= 0.0;
    double r = 1.0;

//    d[0] = 1.0;
//    u[0] = 0.0;
//    b[0] = y0;
    d[0] = (beta-alfa/h);
    u[0] = alfa/h;
    b[0] = -gamma;
    double h_kw=h*h;
    double a=h;
    for(int i=1;i<N-1;i++)
    {

        l[i-1] = p/h_kw-q/(2*h);
        d[i] = -2*p/h_kw+r;
        u[i] = p/h_kw+q/(2*h);
        b[i] = -(2*sin(i*h));
//        l[i-1] = 1.0 / (h_kw);
//        d[i] = (-2.0+h_kw)/h_kw;
//        u[i] = 1.0 / (h_kw);
//        b[i] = -2.0*sin(a);
//        a+=h;

//        l[i-1] = 1.0/h_kw;
//        d[i] = -2/h_kw+1;
//        u[i] = 1.0/h_kw;
//        b[i] = -2*sin((double)i*h);

//        d[i] = -(2.0-h*h);//-h*h);
//        u[i] = 1.0 ;
//        b[i] = -2*sin((double)i*h)*h*h;

//        		 l[i-1] = 1.0 / (h * h);
//		 d[i] = -(8.0 + h * h)/(4.0 * h * h);
//		 u[i] = 1.0 / (h * h);
//		 b[i] = 0.0;
    }
    d[N-1] = fi/h+ni;
    l[N-2] = fi/h;
    b[N-1] = -theta;
//    d[N-1]=1.0;
//    l[N-2]=0.0;
//    b[N-1]=yn;

//    for(int i=0;i<N-1;i++)
//    {
//        printf("%e ",l[i]);
//    }
//    cout << endl<<endl;
//    for(int i=0;i<N;i++)
//    {
//        printf("%e ",d[i]);
//    }
//    cout << endl<<endl;
//    for(int i=0;i<N-1;i++)
//    {
//        printf("%e ",u[i]);
//    }
//    cout << endl<<endl;
//    for(int i=0;i<N;i++)
//    {
//        printf("%e ",b[i]);
//    }
//    cout << endl<<endl;


    thomas(l,u,d,b,x,N);



}
double U(double x)
{
    return x*cos(x);
}

int main()
{
    FILE* plik = fopen("a.txt", "w");
    const int N = 50;
    double x0=0.0;
    double xn=M_PI;
    double h=(xn-x0)/(N-1);
    double *l = new double[N];
    double *d = new double[N];
    double *u = new double[N];
    double *b = new double[N];
    double *x_kon =new double[N];
    dyskretyzacja_konwencjonalna(h,0,-M_PI,N,l,u,d,b,x_kon);
    //cout << x_kon[0]<<endl;
    //cout <<x_kon[N-1]<<endl;

    for(int i=0;i<N;i++)
    {
        //cout << x_kon[i] << "\t" << U(i*h) <<"\t" <<fabs(x_kon[i]-U(i*h))<< endl;
        //printf("%e %e\n",i*h,x_kon[i]);
        fprintf(plik,"%e %e %e\n",i*h,x_kon[i],i*h*cos(i*h));
    }
    fclose(plik);
    return 0;
}

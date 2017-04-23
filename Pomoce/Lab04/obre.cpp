#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const double tolx=0.0000001;
const double tolf=0.0000001;
const int nmax=150;

double jakob[3][3];
double xa[3];
double f[3];

double f1(double x,double y)
{
    return x*x+y*y-1;
}

double f1dx(double x)
{
    return 2*x;
}
double f1dy(double y)
{
    return 2*y;
}
double f1dz(double z)
{
    return 0;
}

double f2(double x,double y)
{
    return 2*x*x-y+(1/2);
}
double f2dx(double x)
{
    return 4*x;
}
double f2dy(double y)
{
    return -1;
}
double f2dz(double z)
{
    return 0;
}

double f3(double x,double y,double z)
{
    return tan(x/4)+y*y*y-z*z*z*z*z;
}
double f3dx(double x)
{
    return tan(x/4)*tan(x/4);
}
double f3dy(double y)
{
    return 3*y*y;
}
double f3dz(double z)
{
    return -5*z*z*z*z;
}

void jacob(double x, double y, double z, double (&jakob)[3][3])
{
    jakob[0][0]=f1dx(x);
    jakob[0][1]=f1dy(y);
    jakob[0][2]=f1dz(z);

    jakob[1][0]=f2dx(x);
    jakob[1][1]=f2dy(y);
    jakob[1][2]=f2dz(z);

    jakob[2][0]=f3dx(x);
    jakob[2][1]=f3dy(y);
    jakob[2][2]=f3dz(z);

}

void xw(double x, double y, double z, double (&xa)[3])
{
    xa[0]=x;
    xa[1]=y;
    xa[2]=z;

}

void fw(double x, double y, double z, double (&f)[3])
{
    f[0]=f1(x,y);
    f[1]=f2(x,y);
    f[2]=f3(x,y,z);
}

double det(double (&jakob)[3][3])
{
    return jakob[0][0]*jakob[1][1]*jakob[2][2]+jakob[0][1]*jakob[1][2]*jakob[2][0]+jakob[0][2]*jakob[1][0]*jakob[2][1]-jakob[0][0]*jakob[1][2]*jakob[2][1]-jakob[0][1]*jakob[1][0]*jakob[2][2]-jakob[0][2]*jakob[1][1]*jakob[2][0];
}

//void dolacz(double )

//void przeksztalc(double (&jakob)[3][3])
//{
    //for(int i=0;i<3;i++)
    //{
       // jakob[1][i]-=2*jakob[2][i];
        //jakob[2][i]-=3*jakob[1][i];
       // jakob[2][i]*=0.5;
   // }
//}

void pomnoz(double x, double (&jakob)[3][3])
{
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            jakob[i][j]*=x;
        }
    }
}

double norma(double a,double b, double c)
{
    return sqrt(a*a+b*b+c*c);
}

void cramer(double (&jakob)[3][3],double &x,double &y,double &z,double a,double b,double c)
{
    double deter=det(jakob);
    double j2[3][3];
    j2[0][0]=a;
    j2[1][0]=b;
    j2[2][0]=c;
    for(int i=0;i<3;i++)
    {
        for(int j=1;j<3;j++)
        {
            j2[i][j]=jakob[i][j];
        }
    }
    double deter2=det(j2);
    j2[0][1]=a;
    j2[1][1]=b;
    j2[2][1]=c;
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            if(j!=1)
                j2[i][j]=jakob[i][j];
        }
    }
    double deter3=det(j2);

    j2[0][2]=a;
    j2[1][2]=b;
    j2[2][2]=c;
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<2;j++)
        {
            j2[i][j]=jakob[i][j];
        }
    }
    double deter4=det(j2);

    x=deter2/deter;
    y=deter3/deter;
    z=deter4/deter;

}


int main()
{
    int n;
    //xw(1,1,1,xa);
    //fw(1,1,1,f);
    //jacob(1,1,1,jakob);
    //przeksztalc(jakob);
    double est=1.0;
    double residuum=1.0;
    double x_,y_,z_;

    xw(1,1,1,xa);
    fw(1,1,1,f);
    jacob(1,1,1,jakob);
    n++;
    printf("n=%d, x=%lf, y=%lf, z=%lf\n",n, xa[0],xa[1],xa[2]);
    while( (n<nmax) && (est>tolf) && (residuum>tolx))
    {


        jacob(xa[0],xa[1],xa[2],jakob);
        cramer(jakob,xa[0],xa[1],xa[2],-f1(xa[0],xa[1]),-f2(xa[0],xa[1]),-f3(xa[0],xa[1],xa[2]));
        //przeksztalc(jakob);
        //xa[n+1]=xa[n]-jakob[3][3]*f[n];
        //xa[0]+=((f[0]-(jakob[0][1]*xa[1]+jakob[0][2]*xa[2]))/jakob[0][0]);
        //xa[1]+=((f[1]-(jakob[1][1]*xa[1]+jakob[1][2]*xa[2]))/jakob[1][0]);
        //xa[2]+=((f[2]-(jakob[2][1]*xa[1]+jakob[2][2]*xa[2]))/jakob[2][0]);
        xa[0]+=x_;
        xa[1]+=y_;
        xa[2]+=z_;
        est=norma(x_,y_,z_);
        residuum=norma(f1(xa[0],xa[1]),f2(xa[0], xa[1]), f3(xa[0], xa[1], xa[2]));
        printf("n=%d, x=%lf, y=%lf, z=%lf, norma=%lf, residuum=%lf\n",n, xa[0],xa[1],xa[2],est,residuum);

        n++;
    }
    //printf("%lf",xa[n+1]);

}

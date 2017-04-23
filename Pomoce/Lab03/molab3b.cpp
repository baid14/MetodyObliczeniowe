    #include <iostream>
    #include <cmath>
    #include <time.h>
    #include <iomanip>
	#include <cstdlib>
	
    using namespace std;
     
    double *metodaThomasa(double **A, double *b,int m)
    {
            double n,r,suma=0;
            double* x = new double [m];
            for(int i=1;i<m;i++)
            {
                    n=A[i-1][i-1];
                    r=b[i-1];
                    A[i][i]=A[i][i]-A[i][i-1]*pow(n,-1)*A[i-1][i];
                    b[i]=b[i]-A[i][i-1]*pow(n,-1)*r; 
                    A[i][i-1]=0;
            }
            cout<<"\nMacierz A po przeksztalceniu:\n"<<endl;                  
            for(int i=0;i<m;i++)
            {
                    for(int j=0;j<m;j++)
                    {
                            if(A[i][j]>=0&&A[i][j]<10) printf("  %.0lf  ",A[i][j]);
                            else printf(" %.0lf  ",A[i][j]);
                    }
                    cout<<endl;
            }
            cout<<"\nwektor b po przeksztalceniu:\n"<<endl;
            for(int j=0;j<m;j++) printf(" %.0lf\n",b[j]);
            
            x[m-1]=b[m-1]/A[m-1][m-1];
            for(int i=m-2;i>=0;i--)
            {
                    for(int j=i+1;j<m;j++)
                    {
                            suma=suma+A[i][j]*x[j];
                    }
                    x[i]=(b[i]-suma)/A[i][i];
                    suma=0;
            }
            return x;
    }
    
    
    int main()
    {
            double n,r,suma=0;
            double b[5] = {18,32,44,36,22},x[5];
            double A[5][5] = {10,4,0,0,0,
                              1,11,3,0,0,
                              0,2,12,2,0,
                              0,0,3,13,1,
                              0,0,0,4,14};
                              
            double* x2 = new double [5];
                              
            double** tab = new double *[5];
            for (int i=0;i<5;i++) tab[i] = new double[5];
            
            double* btab = new double [5];
            
            for(int i=0;i<5;i++) memcpy(tab[i],A[i],sizeof(double)*5);
            memcpy(btab,b,sizeof(double)*5);
            
            //-----------------------------------------------------------------------------------
            cout<<"Macierz A:\n"<<endl;                  
            for(int i=0;i<5;i++)
            {
                    for(int j=0;j<5;j++)
                    {
                            if(A[i][j]>=0&&A[i][j]<10) printf("  %.0lf  ",A[i][j]);
                            else printf(" %.0lf  ",A[i][j]);
                    }
                    cout<<endl;
            }
            cout<<"\nwektor b:\n"<<endl;
            for(int j=0;j<5;j++) printf(" %.0lf\n",b[j]);
            
            cout<<"_______________________________________";
            cout<<endl<<"Rozwiazanie ukladu Ax=b metoda Thomasa:"<<endl;
            
            //-----------------------------------------------------------------------------------
            /*
            for(int i=1;i<5;i++)
            {
                    n=A[i-1][i-1];
                    r=b[i-1];
                    A[i][i]=A[i][i]-A[i][i-1]*pow(n,-1)*A[i-1][i];
                    b[i]=b[i]-A[i][i-1]*pow(n,-1)*r; 
                    A[i][i-1]=0;
            }
            cout<<"\nMacierz A po przeksztalceniu:\n"<<endl;                  
            for(int i=0;i<5;i++)
            {
                    for(int j=0;j<5;j++)
                    {
                            if(A[i][j]>=0&&A[i][j]<10) printf("  %.0lf  ",A[i][j]);
                            else printf(" %.0lf  ",A[i][j]);
                    }
                    cout<<endl;
            }
            cout<<"\nwektor b po przeksztalceniu:\n"<<endl;
            for(int j=0;j<5;j++) printf(" %.0lf\n",b[j]);
            
            cout<<"\nRozwiazanie (wektor x):\n"<<endl;
            x[4]=b[4]/A[4][4];
            for(int i=3;i>=0;i--)
            {
                    for(int j=i+1;j<5;j++)
                    {
                            suma=suma+A[i][j]*x[j];
                    }
                    x[i]=(b[i]-suma)/A[i][i];
                    suma=0;
            }
            for(int j=0;j<5;j++) printf(" %.0lf\n",x[j]);
            cout<<endl;
            */
            //-----------------------------------------------------------------------------------
            
            x2=metodaThomasa(tab,btab,5);
            
            cout<<"\nRozwiazanie (wektor x):\n"<<endl;
            for(int j=0;j<5;j++) printf(" %.3lf\n",x2[j]);
            
            //-----------------------------------------------------------------------------------
            
            cout<<endl;
            		
    system("pause");
    return 0;
    }



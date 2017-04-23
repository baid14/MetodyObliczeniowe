    #include <iostream>
    #include <cmath>
    #include <time.h>
    #include <iomanip>
	#include <cstdlib>
	
    using namespace std;
     
    double **odwrocMacierz(double **tab, int n)
    {
            double temp;
            double *temp2 = new double [n];
            double** tab3 = new double *[n];
            for (int i=0;i<n;i++)
            {
                tab3[i] = new double[n];
                for(int j=0;j<n;j++) tab3[i][j]=0;
            }
            
            double** tab2 = new double *[n];
            for (int i=0;i<n;i++) tab2[i] = new double[2*n];
            
            //do macierzy odwracanej dopisujemy macierz jednostkow¹ I
            for(int i=0;i<n;i++)
            {
                    for(int j=n;j<2*n;j++)
                    {
                            memcpy(tab2[i],tab[i],sizeof(double)*n);
                            if(j%n==i) tab2[i][j]=1;
                            else tab2[i][j]=0;
                    }
            }
            //dzielimy kolejne wiersze przez odpowiadajacy im element przekatnej
            for(int i=0;i<n;i++)
            {
                    if(tab2[i][i]==0)
                    {
                                     for(int k=i;k<n;k++){ 
                                     if(tab2[k][i]!=0){
                                                      temp2=tab2[i];
                                                      tab2[i]=tab2[k];
                                                      tab2[k]=temp2;
                                                      } 
                                     }
                    }
                    
                    temp=tab2[i][i];
                    if(temp!=0) for(int j=0;j<2*n;j++) tab2[i][j]/=temp;
                    else{ 
                         cout<<"zero w "<<i+1<<" wierszu!"<<endl; //return tab3; 
                         }
            
                    //od macierzy odejmujemy odpowiednie elementy
                    //aby wyzerowaæ elementy poza przek¹tnymi macierzy
                    for(int k=0;k<n;k++)
                    {
                            temp=tab2[k][i];
                            if(k!=i) for(int j=0;j<2*n;j++) tab2[k][j]-=(temp*tab2[i][j]);
                    }
            }
            
            for(int i=0;i<n;i++)
            {
                    for(int j=n;j<2*n;j++) tab3[i][j%n]=tab2[i][j];
            }
            return tab3;
    }
    
    
    int main()
    {
            int n;
            cout<<"Podaj rozmiar tablicy: ";
            cin>>n;
            
            double **tabwyn;
            
            double** tab1 = new double *[n];
            for (int i = 0; i < n; ++i)
            {
                    tab1[i] = new double[n];
            }
            
            cout<<"Podaj elementy tablicy: "<<endl;
            for(int i=0;i<n;i++)
            {
                    for(int j=0;j<n;j++) cin>>tab1[i][j];
            }
            
            tabwyn=odwrocMacierz(tab1,n);
            
            cout<<"\nPo odwroceniu:\n"<<endl;
            for(int i=0;i<n;i++)
            {
                    for(int j=0;j<n;j++)
                    {
                            if(tabwyn[i][j]>=0) printf(" %.3lf  ",tabwyn[i][j]);
                            else printf("%.3lf  ",tabwyn[i][j]);
                    }
                    cout<<endl;
            }
            		
    system("pause");
    return 0;
    }



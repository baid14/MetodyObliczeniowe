    #include <cstdlib>
    #include <iostream>
    #include <iomanip>
    #include <cmath>
     
    using namespace std;
     
    /* WYŒWIETLENIE MACIERZY */
    void show(double aa[4][4], int n, int indeks[]){
            for(int i=0; i<n; i++){
                    for(int j=0; j<n; j++){
                            cout << setw(4)<< aa[indeks[i]][j];//setw()- ustawia odstepy pomiedzy kolejnymi wypisywanymi znakami
                    }
            cout << endl;
            }
    }
    /* DEKOMPOZYCJA LU */
    void dek_LU(double aa[4][4], int n, int indeks[]){
            int l, max_w, i, j, k, temp;          
            double max;                  
            for (k=0; k<n; k++){          // sterowanie wielkoœci¹ macierzy do eliminacji Gaussa
                    for(i=k+1; i<n; i++){
                            for(j=n-1; j>=k; j--){
                                    double e=aa[indeks[i]][k]/aa[indeks[k]][k];
                                    if(aa[indeks[k]][k]!=0){
                                            aa[indeks[i]][j]=aa[indeks[i]][j]-aa[indeks[k]][j]*(aa[indeks[i]][k]/aa[indeks[k]][k]);
                                            if(j==k)
                                                    aa[indeks[i]][j]=e;
                                    }
                                    else{
                                            max=aa[indeks[k]][k];
                                            l=k+1;
                                            while(l<n){
                                                    if(fabs(aa[indeks[l]][k])>max){
                                                            max=aa[indeks[l]][k];
                                                            max_w=l;
                                                    }
                                                    l++;
                                            }
                                            if(max_w==n){
                                                    cout << "Brak rozwiazan, lub nieskonczenie wiele rozwiazan.\n";
                                                    system("pause");
                                            }
                                            temp=indeks[k];
                                            indeks[k]=indeks[max_w];
                                            indeks[max_w]=temp;
                                    }
                            }
                    }
            }
    } 
    /* ROZWIAZANIE UKLADU */
    void rozwiazanie(double a[4][4], int n, double b[4], double x[4], double y[4], int indeks[]){
            cout << "Macierz U:\n\n";
            int i, j;
            for(i=0; i<n; i++){
                    for(j=0; j<n; j++){
                            if(j<i)
                                    cout << setw(4) << 0;
                            else
                                    cout << setw(4) << a[indeks[i]][j];
                    }
                    cout << endl;
            }
            cout << endl << endl;
            cout << "Macierz L:\n\n";
            for(i=0; i<n; i++){
                    for(j=0; j<n; j++){
                            if(j>i)
                                    cout << setw(4) << 0;
                            else{
                                    if(j==i)
                                            cout << setw(4) << 1;
                                    else
                                            cout << setw(4) << a[indeks[i]][j];
                            }
                    }
                    cout << endl;
            }
            y[0]=b[indeks[0]]/a[indeks[0]][0];
            double suma=0;
            for(i=1; i<n; i++){
                    for(j=0; j<i; j++)
                            suma=suma+a[indeks[i]][j]* y[j];
                    y[i]=(b[indeks[i]]-suma);
                    suma=0;
            }
            x[n-1]=y[n-1]/a[indeks[n-1]][n-1];
            for(i=n-2; i>=0; i--){
                    for(int j=i+1; j<n; j++)
                            suma=suma+a[indeks[i]][j]*x[j];
                    x[i]=(y[i]-suma)/a[indeks[i]][i];
                    suma=0;
            }
    }
    int main(){
            double A[4][4];
            A[0][0] = 1; A[0][1] = 2; A[0][2] = 2; A[0][3] = 1;    
            A[1][0] = 2; A[1][1] = 4; A[1][2] = 4; A[1][3] = 1;    
            A[2][0] = 2; A[2][1] = 2; A[2][2] = 2; A[2][3] = 1;
            A[3][0] = 1; A[3][1] = 1; A[3][2] = 2; A[3][3] = 1;
            double b[4]={1,2,3,4};
            double y[4], x[4];
            int indeks[4];
            int n = 4, i, j;
            for(i=0; i<n; i++)
                    indeks[i]=i; 
            cout << "Macierz A:\n\n";
            show(A, n, indeks);
            cout << "\nWektor b:\n\n";
            cout << "         b=[ ";
            for(i=0; i<n; i++){
                    cout << b[i] << " ";
            }  
            cout << "]";  
            cout << "\n\n--------------------------------------------------------------------------------\n";  
         
            dek_LU(A, n, indeks);
            rozwiazanie(A,n, b, x, y, indeks);
            cout << "\nRozwiazanie: \n\n" << "x=[ ";
            for(int i=0; i<n; i++){
                             cout << x[i] << "  ";}cout << "]\n\n";
            system("pause");
    }

#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>

/*

Metody Obliczeniowe
By Wahuu

*/


using namespace std;


double **nowaMacierz(int n, int m);
void usunMacierz(double **mac, int n);
void wypelnijMacierz(double **mac, int m, int n);
void wypiszMacierz(double **mac, int m, int n);
int szukaj(double **mac, int n, int j, int *r);
void elmnGauss(double **mac, int *r);
void macierzL(double **l, double *b, int *r, int n);
void macierzU(double **u, double *b, int *p, int n);

int main(int argc, char** argv) {
    double **macierz = nowaMacierz(4, 4);
    double b[4] = {1.0, 2.0, 3.0, 4.0};
    int r[4] = {0, 1, 2, 3};
    int i, j;
    
    cout << "Macierz przed: \n" << endl;
    wypelnijMacierz(macierz, 4, 4);
    wypiszMacierz(macierz, 4, 4);   
    
    elmnGauss(macierz, r);
    
    cout << "Macierz po: \n" << endl;
    wypiszMacierz(macierz, 4, 4);
    
    macierzL(macierz, b, r, 3);
    macierzU(macierz, b, r, 3);
    cout << "Rozwiazanie: \n" << endl;
    
    for(i = 0; i < 4; i++)
        {
        cout << "| " << b[r[i]] << "\t|" <<endl;
        
        }
    
    usunMacierz(macierz, 3); 
    cout<< endl;      
    
    system("pause");
}

/*
 * 
 */

double **nowaMacierz(int n, int m) {     
    double **mac;
    mac = new double *[n];

    for(int i = 0; i < n; i++) 
        mac[i] = new double[m];
    return mac;
}

void usunMacierz(double **mac, int n) {
    for(int i = n - 1; i > 0; i--) 
        delete []mac[i];
    delete []mac;
}

void wypelnijMacierz(double **mac, int m, int n) {
    double tmp[4][4] = { {1.0,2.0,2.0,1.0}, {2.0,4.0,4.0,1.0}, {2.0,2.0,2.0,1.0}, {1.0,1.0,2.0,1.0} };
    for(int i = 0; i < m; i++)
        for(int j = 0; j < n; j++)
            mac[i][j] = tmp[i][j];
}

void wypiszMacierz(double **mac, int m, int n) {
    for(int i = 0; i < m; i++) {
    cout << "| ";
        for(int j = 0; j < n; j++)
            cout << mac[i][j] << "\t";
            cout << "|";
        cout << endl;
    }
}
// szukam rzedu ktorego wartosc bezwzgledna jest najwieksza
int szukaj(double **mac, int n, int j, int *r) {
    int row;    
    for(int i = n + 1; i < j; i++) {
        if(fabs(mac[r[i]][n]) < fabs(mac[r[i+1]][n]))  
            row = r[i+1];
        else
            row = r[i];
    }
    return row;
}


void elmnGauss(double **mac, int *r) {
    int row;                                                     
    double v;    
    for(int k = 0; k < 3; k++) {
        if(mac[r[k]][r[k]] == 0.0) {                                           // sprawdzam czy na pzrekatnej s¹ 0
            row = szukaj(mac, r[k], 3, r);                       
            r[row] = r[k];                                                     // zapisuje w tablicy i podmieniam indeksy
            r[k] = row;                
        }
        
        for(int i = k + 1; i < 4; i++) {
            v = mac[r[i]][k];
            for(int j = k + 1; j < 4; j++)
                mac[r[i]][j] = mac[r[i]][j] - mac[r[k]][j] * (v / mac[r[k]][k]);
            mac[r[i]][k] = v / (mac[r[k]][k]);                                  // Zapisuje L w tej samej macierzy co U
           
        }
    }
}

// rozwiazuje LUx=b      => Ly=b  ->y
void macierzL(double **l, double *b, int *r, int n) {
    double sum = 0;
    for(int i = 0; i <= n; i++) {
        for(int j = 0; j < i; j++)
            sum = sum + l[r[i]][j] * b[r[j]];
        b[r[i]] = b[r[i]] - sum;
        sum = 0;
    }     
}

// rozwiazuje Ux=y -> x
void macierzU(double **u, double *b, int *p, int n) {
     double sum = 0;
     for(int i = n; i >= 0; i--) {
        for(int j = i + 1; j <= n; j++)
           sum = sum + u[p[i]][j] * b[p[j]];
        b[p[i]] = (b[p[i]] - sum) / (u[p[i]][i]);
        sum = 0; 
     }
}

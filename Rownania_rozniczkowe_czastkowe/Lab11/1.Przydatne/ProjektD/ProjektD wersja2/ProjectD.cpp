#include "ProjectD.h"

ProjectD::ProjectD(void)
{
	N = 0;
	M = 0; 
	h = 0;
	dt = 0;
	methName = "";
	D = 1;
	TMAX = 2.0;
	a = 10;
	r = 1.0;
}

ProjectD::~ProjectD(void)
{
}

void ProjectD::info()
{
	cout << endl << "Projekt D";
	cout << endl << "Dyskretyzacja:";
	cout << endl << "   -Klasyczna metoda bezposrednia";
	cout << endl << "   -Metoda posrednia Cranka-Nicolson";
	cout << endl << "Rozwiazanie ukladu rownan:";
	cout << endl << "   -Algorytm Thomasa";
	cout << endl << "Parametry:";
	cout << endl << "d = " << D << "   " << "Tmax = " << TMAX;
	cout << endl << "Przedzial zmiennej przestrzennej: [" << r << ";" << r+a << "]" << endl;
}

void ProjectD::compute(int n)
{
	if ( methName == "CN_T" )
		computeCNT(n);
	else if ( methName == "KMB" )
		computeKMB(n);
	else 
		cout << endl << "Choose method first!";
}

void ProjectD::computeCNT(int n)
{
	cout << endl << "Computing CN + T";
	if ( computed() )
	{
		cout << endl << "Preparing memory...";
		reset();
		cout << "done!";
	}
	N = n;
	h = a / (long double)(N - 1);
	dt = h*h;
	M = static_cast<int>((TMAX / dt) + 1);
	dt = TMAX / (long double)(M-1);
	pos = new long double[N];
	for(int i = 0; i < N; i++)
		pos[i] = r + i*h;
	times = new long double[M];
	for(int i = 0; i < M; i++)
		times[i] = i*dt;

	lambda = D * dt / (h * h);
	lambda2 = lambda / 2.0;
	methName = "CN_T";
	cout << endl << " N = " << N;
    cout << endl << " h = " << h;
    cout << endl << " M = " << M;
    cout << endl << " dt = " << dt;
    cout << endl << " lambda = " << lambda;
	cout << endl << "Start...";
	CElapsed timer;
	timer.Begin();

	x = new long double*[M];
	for (int i = 0; i < M; i++)
		x[i] = new long double[N];
	
	long double *b;
	b = new long double[N];			//wektor wyrazow wolnych w ukladzie rownan Ax = b;

	long double *l;						   //wektor zawierajacy wartosci podprzekatnej macierzy A
	l = new long double[N-1];
	long double *d;						   //wektor zawierajacy wartosci przekatnej macierzy A
	d = new long double[N];
	long double *u;						   //wektor zawierajacy wartosci nadprzekatnej macierzy A
	u = new long double[N-1];
	long double *ptrX;

	ptrX = x[0];
	long double modR = r/h;

	initialValues(ptrX);
	//initialValues(b);

	/*for(int i = 0; i < N; i++)			//rozwiazanie dla poziomu czasowego t=0 to warunek poczatkowy
		ptrX[i] = -b[i];*/

	diskretCNT(l, d, u);				//ladujemy do wektorow odpowiednie wartosci
	
	thomasReduce(l, d, u);				//pierwszy etap algorytmu Thomasa - zerujemy podprzekatna

	for( int j = 1; j < M; j++ ) 
	{ 
		for(int i = 1; i < N-1; i ++)
			b[i] = (-lambda2 + lambda2/(i+modR))* ptrX[i-1] + (lambda - 1.0)*ptrX[i] - (lambda2 + lambda2/(i+modR))*ptrX[i+1];
		b[0] = 0.0;						//warunek brzegowy dla poczatku przedzialu
		b[N-1] = 1.0 - (r/(r+a) * erfc(a/(2.0*sqrt(D*times[j]))));					//warunek brzegowy dla konca przedzialu

		ptrX = x[j];
		thomasSolve(l, d, u, ptrX, b);	//drugi etap Thomasa - wyliczamy rozwiazania na kolejnych 
										//poziomach czasowych (iteracja po j)
										//za kazdym razem aktualizujemy wektor b jak ponizej
	}


	delete [] b;

	delete [] l;
	delete [] d;
	delete [] u;

	elapsed = timer.End();
	cout << "done: " << elapsed << "s" << endl;
}

void ProjectD::computeKMB(int n)
{
	cout << endl << "Computing KMB";
	if ( computed() )
	{
		cout << endl << "Preparing memory...";
		reset();
		cout << "done!";
	}
	N = n;
	h = a / (long double)(N - 1);
	dt = 0.4*h*h;
	M = static_cast<int>((TMAX / dt) + 1);
	dt = TMAX / (long double)(M-1);
	pos = new long double[N];
	for(int i = 0; i < N; i++)
		pos[i] = r + i*h;
	times = new long double[M];
	for(int i = 0; i < M; i++)
		times[i] = i*dt;

	lambda = D * dt / (h * h);
	methName = "KMB";
	cout << endl << " N = " << N;
    cout << endl << " h = " << h;
    cout << endl << " M = " << M;
    cout << endl << " dt = " << dt;
    cout << endl << " lambda = " << lambda;
	cout << endl << "Start...";
	CElapsed timer;
	timer.Begin();
	long double modR = r/h;

	x = new long double*[M];
	for (int i = 0; i < M; i++)
		x[i] = new long double[N];

	long double *ptrX;

	ptrX = x[0];
	initialValues(ptrX);

	for( int j = 1; j < M; j++ ) 
	{
		x[j][0] = 0.0;
		x[j][N-1] = 1.0- (r/(r+a) * erfc(a/(2.0*sqrt(times[j]))));
		for(int i = 1; i < N-1; i++)
			x[j][i] = (lambda-lambda/(i+modR))*ptrX[i-1] + (1.0-2.0*lambda)*ptrX[i] + (lambda+lambda/(i+modR))*ptrX[i+1];
				
		ptrX = x[j];
	}

	elapsed = timer.End();
	cout << "done: " << elapsed << "s" << endl;
}

void ProjectD::setMethod(int meth)
{
	if ( meth == 1 )
		methName = "CN_T";
	else if ( meth == 2 )
		methName = "KMB";
}

void ProjectD::initialValues(long double *b)
{
	for( int i = 0; i < N; i++)
			b[i] = 1.0;
}

void ProjectD::diskretCNT(long double *l, long double *d, long double *u)
{
	long double modR = r/h;
	for( int i = 0; i < N-2; i++ )      // wektor l
         l[i] = lambda2 - lambda2/(i+modR);
    l[N-2]=0.0;

    for( int i = 1; i < N-1; i++ )      // wektor d
         d[i] = -(1.0 + lambda);
    d[N-1] = 1.0;
    d[0] = 1.0;

    u[0]=0.0;
	for( int i = 1; i < N-1; i++ )        // wektor u
         u[i] = lambda2 + lambda2/(i-1+modR); 
}

void ProjectD::thomasReduce(long double *l, long double *d, long double *u)
{
	for( int i = 1; i < N; i++)
		d[i] -= ( l[i-1] / d[i-1] ) * u[i-1];
}
void ProjectD::thomasSolve(long double *l, long double *d, long double *u, long double *x, long double *b)
{
	for( int i = 1; i < N; i++)
		b[i] -= ( l[i-1] / d[i-1] ) * b[i-1];

	x[N-1] = b[N-1] / d[N-1];

	for( int i = N-2; i >= 0; i--)
		x[i] = ( b[i] - u[i] * x[i+1] ) / d[i];
}

void ProjectD::reset()
{
	for(int i = 0; i < M; i++)
		delete [] x[i];
	delete []x ;
	N = 0;
	M = 0; 
	h = 0;
	dt = 0;
	methName = "";
}
#define OUTPUT_DIR "resultsss/"

#include "stdafx.h"
#include <sstream>
#include <direct.h>

void graph1(ProjectD proj, int method);
void graph2(ProjectD proj, long double time, int freq);
void graph3(ProjectD proj, int freq);

int _tmain(int argc, _TCHAR* argv[])
{
	ProjectD proj;
	proj.info();
	
	/*
	 * Obliczenia dla CN + T
	 */
	/*graph1(proj, 1);
	proj.setMethod(1);
	proj.compute(200);
	graph2(proj, 0.05, 20);
	graph2(proj, 0.35, 20);
	graph2(proj, 0.45, 20);
	graph2(proj, 0.73, 20);
	graph2(proj, 0.9, 20);
	graph3(proj, 1200);*/

	/*
	 * Obliczenia dla KMB
	 */
	graph1(proj, 2);
	proj.setMethod(2);
	proj.compute(200);
	graph2(proj, 0.05, 20);
	graph2(proj, 0.35, 20);
	graph2(proj, 0.45, 20);
	graph2(proj, 0.73, 20);
	graph2(proj, 0.9, 20);
	graph3(proj, 1200);

	system("PAUSE");
	return 0;
}

void graph1(ProjectD proj, int method)
{
	cout << endl << "Generating graph 1";
	int maxTests = 1600;
	proj.setMethod(method);
	std::string fileName = OUTPUT_DIR + proj.getName() + "-graph1.txt";

	_mkdir(OUTPUT_DIR);
	ofstream fGraph(fileName.c_str());
	if ( fGraph.is_open() )
	{
		long double **x;
		long double *ptrX;
		int N;
		int M;
		long double h;
		long double *times;

		long double time;
		long double err;
		long double max;

		long double h1 = 0.0;
		long double h2 = 0.0;
		long double Fh1 = 0.0;
		long double Fh2 = 0.0;

		for(int k = 100; k <= maxTests; k*= 2)
		{
			proj.compute(k);
			N = proj.getN();
			h = proj.getH();
			M = proj.getM();
			x = proj.getX();
			times = proj.getTimes();

			ptrX = x[M-1];
			max = 0.0;
			time = times[M-1];
			for(int i = 0; i < N; i++)
			{
				err = ptrX[i] - proj.exactSolution(M-1, i);
				err = fabs( err );
				if ( err > max )
					max = err;
			}
			fGraph << h << "\t" << setiosflags(ios::fixed) << std::setprecision(7) << max << "\t";
			fGraph << proj.getElapsed() << endl;
			if ( k == 200 )
			{
				h1 = h;
				Fh1 = max;
			}
			if ( k == 400 )
			{
				h2 = h;
				Fh2 = max;
			}
		}

		long double dokl;
		dokl = (log(Fh2) - log(Fh1)) / (log(h2) - log(h1));
		cout << endl << "Dokladnosc: " << dokl;
		fGraph << endl << dokl;
		fGraph.close();
		cout << endl << "Graph1 done!";
	}
	else
		cout << endl << "Error opening file for write: " << fileName << "!";
}

void graph2(ProjectD proj, long double time, int freq)
{
	if ( proj.computed() )
	{
		if ( time >= 0 && time <= 1)
		{
			int j = static_cast<int>(time * proj.getM());
			long double dt = proj.getDt();
			time = j*dt;
			int N = proj.getN();
			long double h = proj.getH();
			long double a = proj.getA();
			long double r = proj.getR();
			long double *pos;
			pos = proj.getPos();
			long double *times;
			times = proj.getTimes();
			cout << endl << "Generating graph 2";
			cout << endl << " time = " << time;
			cout << endl << " method: " << proj.getName();
			cout << endl << "Start...";
			std::ostringstream ostime;
			ostime << time;
			string fileName;
			fileName += OUTPUT_DIR;
			fileName += proj.getName();
			fileName += "-" + ostime.str();
			std::ostringstream osfreq;
			osfreq << freq;
			fileName += "-" + osfreq.str();
			fileName += "-graph2.txt";

			long double **x;
			x = proj.getX();
			long double *ptrX;
			ptrX = x[j];

			_mkdir(OUTPUT_DIR);
			ofstream fGraph(fileName.c_str());
			if ( fGraph.is_open() )
			{
				for(int i = 0; i < N; i+= freq)
					fGraph << pos[i] << "\t" << setiosflags(ios::fixed) << std::setprecision(7) << ptrX[i] << endl; 
				fGraph.close();
				cout << "done!";
			}
			else
				cout << endl << "error opening file for write: " << fileName;
		}
		else
			cout << endl << "Usage: 0 <= timer <= 1!";
	}
	else
		cout << endl << "Compute first!";
	cout << endl;
}

void graph3(ProjectD proj, int freq)
{
	if ( proj.computed() )
	{
		cout << endl << "Generating graph 3";
		cout << endl << " method: " << proj.getName();
		cout << endl << "Start...";
		long double dt = proj.getDt();
		int N = proj.getN();
		int M = proj.getM();
		long double *pos;
		pos = proj.getPos();
		long double *times;
		times = proj.getTimes();
		string fileName;
		fileName += OUTPUT_DIR;
		fileName += proj.getName();
		std::ostringstream osfreq;
		osfreq << freq;
		fileName += "-" + osfreq.str();
		fileName += "-graph3.txt";

		_mkdir(OUTPUT_DIR);
		ofstream fGraph(fileName.c_str());
		if ( fGraph.is_open() )
		{
			long double **x;
			x = proj.getX();
			long double *ptrX;

			long double err;
			long double max;

			for( int j = 0; j < M; j+= freq)
			{
				max = 0;
				ptrX = x[j];
				for(int i = 0; i < N; i++)
				{
					err = ptrX[i] - proj.exactSolution(j+1, i);
					err = fabs( err );
					if ( err > max )
						max = err;
				} 
				fGraph << times[j] << "\t" << setiosflags(ios::fixed) << std::setprecision(7) <<  max << endl; 
			}
			fGraph.close();
			cout << "done!";
		}
		else
			cout << endl << "error opening file for write: " << fileName;
	}
	else
		cout << endl << "Compute first!";
	cout << endl;
}
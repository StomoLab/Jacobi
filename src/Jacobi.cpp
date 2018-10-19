//============================================================================
// Name        : Jacobi.cpp
// Author      : T. Suzuki
// Date        : 2018/10/18
//============================================================================

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cmath>

using namespace std;

int main(const int argc, const char *argv[])
{
	if (argc < 2)
	{
		cerr << "Usage: a.out [M]" << endl;
		exit(EXIT_FAILURE);
	}

	const unsigned long m = atoi(argv[1]);  // size of matrix

	#ifdef DEBUG
	cout << "m = " << m << endl << endl;
	#endif

	double *a = new double [m*m];
	double *ad = new double [m*m];

	// Create random symmetric matrix
	srand(time(NULL));
	for (unsigned long i=0; i<m; i++)
		for (unsigned long j=0; j<=i; j++)
			a[ i + m*j ] = a[ j + m*i ] = ad[i+m*j] = ad[j+m*i] = (double)rand() / RAND_MAX;


	#ifdef DEBUG

	for (unsigned long i=0; i<m; i++)
	{
		for (unsigned long j=0; j<m; j++)
			cout << a[ i + m*j ] << ", ";
		cout << endl;
	}
	cout << endl;
	#endif

	const unsigned int p=3;
	const unsigned int q=6;
	assert((p<=m) && (q<=m));
	assert(p<=q);

	double th = M_PI/4.0;
	double c = 1.0;
	double s = 0.0;
	double c2 = 1.0;
	double s2 = 0.0;

	if (a[p+m*p] != a[q+m*q])
	{
		th = atan( -2.0*a[p+m*q] / ( a[p+m*p] - a[q+m*q] ) );

		c2 = cos(th);
		s2 = sin(th);

		th /= 2.0;
		c = cos(th);
		s = sin(th);
	}
	cout << "theta =" << th << endl;

	double tmp = (a[p+m*p] - a[q+m*q])*s2/2.0 + a[p+m*q]*c2;
//	cout << "ad[" << p << "," << q << "] = ad["  << q << "," << p << "] =" << tmp << endl << endl;

	for (unsigned long i=0; i<m; i++)
	{
		if (i==p)
		{
			for (unsigned long j=0; j<m; j++)
			{
				if (j==p) // a[p][p]
				{
					ad[i+m*j] = a[p+m*p]*c*c + a[q+m*q]*s*s - 2.0*a[p+m*q]*s*c;
					#ifdef DEBUG
					cout << "ad[" << i << "," << j << "]\n";
					#endif
				}
				else if (j==q) // a[p][q]
				{
					ad[i+m*j] = tmp;
					#ifdef DEBUG
					cout << "ad[" << i << "," << j << "] = " << tmp << "\n";
					#endif
				}
				else // a[p][j], a[j][p]
				{
					ad[i+m*j] = ad[j+m*i] = a[p+m*j]*c - a[q+m*j]*s;
					#ifdef DEBUG
					cout << "ad[" << i << "," << j << "], ad[" << j << "," << i << "]\n";
					#endif
				}
			}
		}
		else if (i==q)
		{
			for (unsigned long j=0; j<m; j++)
			{
				if (j==p) // a[q][p]
				{
					ad[i+m*j] = tmp;
					#ifdef DEBUG
					cout << "ad[" << i << "," << j << "] = " << tmp << "\n";
					#endif
				}
				else if (j==q) // a[q][q]
				{
					ad[i+m*j] = a[p+m*p]*s*s + a[q+m*q]*c*c + 2.0*a[p+m*q]*s*c;
					#ifdef DEBUG
					cout << "ad[" << i << "," << j << "]\n";
					#endif
				}
				else // a[q][j], a[j][q]
				{
					ad[i+m*j] = ad[j+m*i] = a[p+m*j]*s + a[q+m*j]*c;
					#ifdef DEBUG
					cout << "ad[" << i << "," << j << "], ad[" << j << "," << i << "]\n";
					#endif
				}
			}
		}
	}

	#ifdef DEBUG
	for (unsigned long i=0; i<m; i++)
	{
		for (unsigned long j=0; j<m; j++)
			cout << ad[i+m*j] << ", ";
		cout << endl;
	}
	cout << endl;
	for (unsigned long i=0; i<m; i++)
	{
		for (unsigned long j=0; j<m; j++)
			cout << a[i+m*j] - ad[i+m*j] << ", ";
		cout << endl;
	}
	cout << endl;
	#endif

	delete[] a;
	delete[] ad;

	return EXIT_SUCCESS;
}

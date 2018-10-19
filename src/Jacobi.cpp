//============================================================================
// Name        : Jacobi.cpp
// Author      : T. Suzuki
// Date        : 2018/10/18
//============================================================================

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <cassert>

using namespace std;

/**
 * Frobenius norm
 *   L1 BLAS の dnrm2() を使うべき
 */
double Fnorm( unsigned long m, double* a )
{
	double tmp = 0.0;

	for(unsigned long i=0; i<m*m; i++)
		tmp += a[i]*a[i];
	return sqrt(tmp);
}

int main(const int argc, const char *argv[])
{
	if (argc < 4)
	{
		cerr << "Usage: a.out M p q" << endl;
		exit(EXIT_FAILURE);
	}

	const unsigned long m = atoi(argv[1]);  // size of matrix
	const unsigned long p = atoi(argv[2]);  // a[p][q], a[q][p] are bunished by Givens rotation
	const unsigned long q = atoi(argv[3]);
	assert((p<=m) && (q<=m));
	assert(p<=q);

	// Display size of matrix, element index
	#ifdef DEBUG
	cout << "m = " << m << endl;
	cout << "(p,q) = (" << p << "," << q << ")\n\n";
	#endif

	double *a = new double [m*m];   // Original matrix
	double *ad = new double [m*m];  // ad = g^{T} a g
	double *g = new double [m*m];   // Transformation matrix

	// Create random symmetric matrix
	srand(time(NULL));
	for (unsigned long i=0; i<m; i++)
		for (unsigned long j=0; j<=i; j++)
			a[ i + m*j ] = a[ j + m*i ] = ad[i+m*j] = ad[j+m*i] = (double)rand() / RAND_MAX;

	// Display matrix elements
	#ifdef DEBUG
	cout << "Matrix a:\n";
	for (unsigned long i=0; i<m; i++)
	{
		for (unsigned long j=0; j<m; j++)
			cout << a[ i + m*j ] << ", ";
		cout << endl;
	}
	cout << endl;
	#endif

	double th;      // theta
	double c, s;    // cos(theta), sin(theta)
	double c2, s2;  // cos(2*theta), sin(2*theta)

	if (a[p+m*p] != a[q+m*q])
		th = atan( -2.0*a[p+m*q] / ( a[p+m*p] - a[q+m*q] ) ) / 2.0;
	else
		th = M_PI/4.0;

	#ifdef DEBUG
	cout << "theta =" << th << endl;
	#endif

	c = cos(th);  s = sin(th);
	c2 = cos(2.0*th);  s2 = sin(2.0*th);

	double tmp = (a[p+m*p] - a[q+m*q])*s2/2.0 + a[p+m*q]*c2;

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

	// Display matrix emelemts
	#ifdef DEBUG
	cout << "\nMatrix ad:\n";
	for (unsigned long i=0; i<m; i++)
	{
		for (unsigned long j=0; j<m; j++)
			cout << ad[i+m*j] << ", ";
		cout << endl;
	}
	cout << "\nMatrix a - ad:\n";
	for (unsigned long i=0; i<m; i++)
	{
		for (unsigned long j=0; j<m; j++)
			cout << a[i+m*j] - ad[i+m*j] << ", ";
		cout << endl;
	}
	cout << endl;
	#endif

	cout << "Norm(a) = " << Fnorm(m,a) << endl;
	cout << "Norm(ad) = " << Fnorm(m,ad) << endl;

	// Check routines
	#ifdef DEBUG
	for (unsigned long i=0; i<m; i++)
		for (unsigned long j=0; j<m; j++)
			g[i+m*j] = (i==j) ? 1.0 : 0.0;
	g[p+m*p] = g[q+m*q] = c;
	g[p+m*q] = s; g[q*m*p] = -s;

	#endif

	delete[] a;
	delete[] ad;
	delete[] g;

	return EXIT_SUCCESS;
}

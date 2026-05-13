#include <iostream>
#include <cmath>
#include <vector>
#include <functional>

// Factorial function
double factorial(int x)
{
    if(x == 0 || x == 1) return 1;
    return x * factorial(x - 1);
}

// Double factorial
double double_factorial(int n)
{
	if(n == 0 || n == 1)
		return 1;

	return double_factorial(n - 2)*n;
}

// Power with integer exponent
double intpower(double v, int k)
{
	double val = 1.0;
	
	for(int i = 1; i <= k; ++i)
		val = val*v;
	
	return val;	
}

// Flattening and mapping a 2D grid onto a 1D vector linearly, row by row.
int idx(int n, int l, int lmax, int nmin) 
{
    return (n - nmin) * (lmax + 1) + l; 
}

// Legendre polynomials
double legendre_P(int order, double x)
{
	// P_0 (x) = 1
    if (order == 0)
        return 1.0;

	// P_1 (x) = x
    if (order == 1)
        return x;

	// P_n (x) for n ≥ 2
	// Using the recurrence relation: P_n(x) = [(2n - 1) x P_{n-1}(x) - (n - 1) P_{n-2}(x)] / n
	double P_n;			  // P_{n} (x)
	double Pminus2 = 1.0; // P_{n-2} (x) -- starts as P_{0} (x) = 1.0
    double Pminus1 = x;   // P_{n-1} (x) -- starts as P_{1} (x) = x

    for (int n = 2; n <= order; ++n)
    {
        P_n = ((2.0 * n - 1.0) * x * Pminus1 - (n - 1.0) * Pminus2) / n; // current P_n (x)

        // Updates for next loop
		Pminus2 = Pminus1; // P_{n-2} -> P_{n-1}
        Pminus1 = P_n;     // P_{n-1} -> P_{n}
    }

    return P_n;
}
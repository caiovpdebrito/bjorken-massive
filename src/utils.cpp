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

#pragma once
#include <vector>
#include <cmath>
#include <iostream>

// Factorial
double factorial(int x);

// Double factorial
double double_factorial(int n);

// Power with integer exponent 
double intpower(double v, int k);

// Index mapping function
int idx(int n, int l, int lmax, int nmin);

// Legendre polynomials 
double legendre_P(int order, double x);
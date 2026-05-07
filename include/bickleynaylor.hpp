#ifndef BICKLEYNAYLOR_HPP
#define BICKLEYNAYLOR_HPP

#include "constants.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

using constants::pi;

// (Exponential) mapping function
double phi(double t);

// Derivative of the mapping function
double dphi(double t); 

// Integrand of the Bickley-Naylor function
double integrand(int n, double x, double w);

// Computation of the Bickley-Naylor function using the Tanh-Sinh quadrature
double Ki_tanh_sinh(int n, double x,
                    double tol = 1e-14,
                    double h = 0.1,
                    int maxK = 60);

// Generate output data file with ĸ_n (z)
void write_Ki_table(const std::string& filename,
                 int nMin,
                 int nMax,
                 double zMin,
                 double zMax,
                 int numPoints);

#endif
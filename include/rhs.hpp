#ifndef RHS_HPP
#define RHS_HPP

#include <vector>
#include <cmath>
#include <iostream>

class KiTable;

// Moments of the equilibrium single-particle distribution function
double equilibrium_moments(int n, 
                           double mass, 
                           double temperature, 
                           double alpha,
                           const KiTable& BickleyNaylor); // ĸ_n(z))

// RHS function
void rhs(double tau,
         double eta_over_s, // shear viscosity rescaled by the entropy density 
         double mass,
         const std::vector<double>& rho,
         std::vector<double>& drhodtau,
         int nmin,
         int nmax,
         int lmax,
         const KiTable& BickleyNaylor);

#endif
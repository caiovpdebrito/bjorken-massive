#include "rhs.hpp"
#include "utils.hpp"
#include "KiTable.hpp"
#include "bickleynaylor.hpp"
#include "constants.hpp"

#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <fstream>
#include <stdexcept>
#include <string>

using constants::pi;

//==========================================//
//======= COMPUTE MOMENTS OF f_{0k} ========//
//==========================================//

double equilibrium_moments(int n, 
                           double mass, 
                           double T, // temperature
                           double alpha, // thermal potential
                           const KiTable& BickleyNaylor) // ĸ_n(z)
{
    double z = mass/T;

    // Massless limit
    if(z < 1e-10)
        return std::exp(alpha)/(2*pi*pi)
               * std::pow(T,n+2)
               * factorial(n+2);

    // General expression 
    return std::exp(alpha)/(2*pi*pi)
           * std::pow(T,n+2)
           * BickleyNaylor.kappa(n, z); 
}

//==========================================//
//==== DEFINE COEFFICIENTS IN THE ODEs =====//
//==========================================//

// General contributions 
double A_coeff(int n, int l) { return -2*l*(2*l-1)*(n+2*l) / ((4*l+1)*(4*l-1)); }
double B_coeff(int n, int l) { return -((2*l*(2*l+1) + n*(24*l*l+12*l-3)) / (3*(4*l+3)*(4*l-1)) - 2/3); } // { return -(2*l*(2*l+1)*(2*n+3) - (n+2)) / ((4*l+3)*(4*l-1)); }
double C_coeff(int n, int l) { return -(n-2*l-1)*(2*l+2)*(2*l+1) / ( (4*l+1)*(4*l+3) ); }

// Intrinsically massive constributions 
double D_coeff(int n, int l) { return 4*l*(n-2*l-1)*((4*l+1) - 1) / ( 3*(4*l-1)*(4*l+3) ); }
double E_coeff(int n, int l) { return -2*l*(2*l-1)*(n-2*l-1) / ( (4*l-1)*(4*l+1) ); }
double F_coeff(int n, int l) { return 2*l*(2*l-1)*(2*n-1) / ( (4*l-1)*(4*l+1) ); }

//==========================================//
//======== COMPUTE RHS OF THE ODEs =========//
//==========================================//

void rhs(double tau,
         double eta_over_s, // shear viscosity rescaled by the entropy density
         double mass, // particle mass in [GeV]
         const std::vector<double>& rho,
         std::vector<double>& drhodtau,
         int nmin,
         int nmax,
         int lmax,
         const KiTable& BickleyNaylor)
{
    // Logarithm time for better numerical behavior
    double t = std::exp(tau); 

    int N = (nmax - nmin + 1) * (lmax + 1); // Allocate the size of the vector
    drhodtau.assign(N, 0.0); 

    // Compute thermodynamic quantities within Landau matching conditions
    double n0 = rho[idx(1, 0, lmax, nmin)]; // particle density
    double e0 = rho[idx(2, 0, lmax, nmin)]; // energy density
    double temperature = BickleyNaylor.T_Landau(mass, e0/n0); 
    double alpha = BickleyNaylor.alpha_Landau(mass/temperature, temperature, n0); 

    // Compute relaxation time -- assume \tau_R = 5\eta_s/(\varepsilon_0 + P_0)
    double tau_R = 5 * eta_over_s * (1 - alpha/4)/temperature; // shear relaxation time
    
    for(int i = nmin; i <= nmax; ++i){
        for(int j = 0; j <= lmax; ++j){
            int k = idx(i, j, lmax, nmin); // index of the flat vector 
            double val = 0.0;

            // coefficients for the pair of indices (i,j)
            double A = A_coeff(i, j);
            double B = B_coeff(i, j);
            double C = C_coeff(i, j);
            double D = D_coeff(i, j);
            double E = E_coeff(i, j);
            double F = F_coeff(i, j);

            // Couplings in rank (rho_{n, l-1}, rho_{n, l}, rho_{n, l+1}) 
            if(j > 0) val += A * rho[idx(i, j-1, lmax, nmin)]; // exclude moments of negative ranks
            val += B * rho[k];
            if(j < lmax) val += C * rho[idx(i, j+1, lmax, nmin)]; // do no account for moments outside truncation in rank

            // Couplings in energy -- exclude terms outside of the truncation
            if(i - 2 >= nmin) val += intpower(mass, 2) * D * rho[idx(i-2, j, lmax, nmin)];
            if(i - 2 >= nmin && j > 0) val += intpower(mass, 4) * E * rho[idx(i-2, j-1, lmax, nmin)];
            if(i - 4 >= nmin && j > 0) val += intpower(mass, 2) * F * rho[idx(i-4, j-1, lmax, nmin)];

            // Collision term (RTA) -- the equilibrium moments for n=1 and n=2 are zero (Landau matching conditions)
            if(i !=1 && i!= 2 && j == 0) val -= (t / tau_R) * (rho[k] - equilibrium_moments(i, mass, temperature, alpha, BickleyNaylor));

            drhodtau[k] = val;
        }
    }
}

#include "bickleynaylor.hpp"
#include "constants.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

using constants::pi;

// ====================================================== //
//              Mapping function w = phi(t)               //
// ====================================================== //

double phi(double t) { return std::exp(pi/2 * std::sinh(t)); }

// ====================================================== //
//          Derivative of the mapping dphi(t)/dt          //
// ====================================================== //

double dphi(double t) { return phi(t) * pi/2 * std::cosh(t); }

// =========================================== //
//  Integrand of the Bickley-Naylor functions  //
// =========================================== //

double integrand(int n, double x, double w)
{
    double exponent = -x * std::cosh(w);

    // Underflow guard
    if (exponent < -745.0)
        return 0.0;

    return std::exp(exponent) * std::pow(std::cosh(w), -n); // Bickley-Naylor functions
    // return std::exp(exponent) * std::pow(std::cosh(w), n) * std::pow(std::sinh(w), 2);
}

// ====================================================== //
//          Compute Ki using Tanh-Sinh quadrature         //
// ====================================================== //

double Ki_tanh_sinh(int n, double x,
                    double tol,
                    double h,
                    int maxK)
{
    double sum = 0.0;

    for (int k = -maxK; k <= maxK; ++k) // symmetric sum range for simplicity
    {
        
        double t = k * h;
        double contrib = integrand(n, x, phi(t)) * dphi(t);

        sum += contrib;

        double weighted = h * contrib;

        // Adaptive truncation: stop once contribution to integral is negligible
        if (k > 8 && std::abs(weighted) < tol)
            break;
    }

    return h * sum;
}

// ====================================================== //
//            Write output data file with ĸ_n             //
// ====================================================== //

void write_Ki_table(const std::string& filename,
                 int nMin,
                 int nMax,
                 double xMin,
                 double xMax,
                 int numPoints)
{
    std::ofstream file(filename);
    file << std::scientific << std::setprecision(12);

    file << "x";
    for (int n = nMin; n <= nMax; ++n)
        file << "\t" << n;
    file << "\n";

    for (int i = 0; i <= numPoints; ++i)
    {
        double x = xMin + (xMax - xMin) * i / numPoints;

        file << x;

        for (int n = nMin; n <= nMax; ++n)
        {
            double kappa_n =
                std::pow(x, n + 2)
                * (Ki_tanh_sinh(-n - 2, x) - Ki_tanh_sinh(-n, x));

            file << "\t" << kappa_n;
        }

        file << "\n";
    }
}


// ###################################################### //
// ###################################################### //

// int main()
// {

//     // Truncation
//     const int nMin = -10;
//     const int nMax = 10;

//     // Numerical limits
//     const double zMin = 1e-3;
//     const double zMax = 1e3;
//     const int numPoints = 1000;

//     std::string filename = "data/Ki_table.dat";

//     write_Ki_table(filename,
//                 nMin,   // nMin
//                 nMax,   // nMax
//                 zMin,
//                 zMax,
//                 numPoints);

//     std::cout << "Table written to Ki_table.dat\n";

//     return 0;
// }

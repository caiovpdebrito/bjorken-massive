#include "rk4.hpp"
#include "rhs.hpp"
#include "KiTable.hpp"
#include "utils.hpp"
#include "constants.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using constants::GeVfm;

double relaxation_time(double eta_over_s,
                       double e0,
                       double P0,
                       double s0)
{
    if (e0 + P0 <= 0.0)
        throw std::runtime_error("ERROR: e0 + P0 must be positive.");

    if (s0 <= 0.0)
        throw std::runtime_error("ERROR: entropy density became non-positive.");

    return 5.0 * eta_over_s * s0 / (e0 + P0);
}

double scaled_time(double tau,
                   double eta_over_s,
                   double mass,
                   const std::vector<double>& rho,
                   int nmin,
                   int lmax,
                   const KiTable& ki_table)
{
    double t = std::exp(tau);

    double n0 = rho[idx(1, 0, lmax, nmin)];
    double e0 = rho[idx(2, 0, lmax, nmin)];

    double T = ki_table.T_Landau(mass, e0 / n0);
    double alpha = ki_table.alpha_Landau(mass / T, T, n0);

    double P0 = n0 * T; // classical ideal gas pressure
    double s0 = (e0 + P0 - alpha * T * n0) / T;

    double tau_R = relaxation_time(eta_over_s, e0, P0, s0);

    return t / tau_R;
}

void write_observables(std::ofstream& file,
                       double tau,
                       double eta_over_s,
                       double mass,
                       const std::vector<double>& rho,
                       int l,
                       int nmin,
                       int nmax,
                       int lmax,
                       const KiTable& ki_table)
{
    double t = std::exp(tau);

    double n0 = rho[idx(1, 0, lmax, nmin)];
    double e0 = rho[idx(2, 0, lmax, nmin)];

    double T = ki_table.T_Landau(mass, e0 / n0);
    double alpha = ki_table.alpha_Landau(mass / T, T, n0);

    double P0 = n0 * T; // classical ideal gas pressure
    double s0 = (e0 + P0 - alpha * T * n0) / T;

    double tau_R = relaxation_time(eta_over_s, e0, P0, s0);

    file << t/tau_R;

    for (int n = nmin; n <= nmax; ++n)
    {
        double rho_nl = rho[idx(n, l, lmax, nmin)];
        double rho_eq_n0 = equilibrium_moments(n, mass, T, alpha, ki_table);

        file << "\t" << std::abs(rho_nl) / rho_eq_n0;
    }

    file << "\n";
}

int main(int argc, char* argv[])
{
    if (argc != 5)
    {
        std::cerr << "Usage: " << argv[0] << " <l> <w0> <wf> <nsteps>\n";
        return 1;
    }

    // Gather arguments parsed 
    int l_plot = std::stoi(argv[1]);
    double w0  = std::stod(argv[2]);
    double wf  = std::stod(argv[3]);
    int nsteps = std::stoi(argv[4]);

    //==========================================//
    //============= PHYSICAL INPUT =============//
    //==========================================//

    // Truncation
    int nMin = 0;
    int nMax = 40;
    int lMax = 10;

    int N = (nMax - nMin + 1) * (lMax + 1); // define dimension of the rho vector
    std::vector<double> rho(N, 0.0); // allocate vector of corresponding size

    // System properties
    double eta_over_s = 1.0;
    double mass = 0.135 * GeVfm;

    // Initial conditions
    double initialTemp = 0.2 * GeVfm;
    double alpha0 = 0.0;

    //==========================================//
    //================= GUARDS =================//
    //==========================================//

    // Check if chosen l is within the truncation
    if (l_plot < 0 || l_plot > lMax)
    {
        std::cerr << "ERROR: l must satisfy 0 <= l <= " << lMax << "\n";
        return 1;
    }

    // Check if initial and final times are positive
    if (w0 <= 0.0 || wf <= 0.0)
    {
        std::cerr << "ERROR: w0 and wf must be positive.\n";
        return 1;
    }

    // Safety: final time is smaller than initial time
    if (wf <= w0)
    {
        std::cerr << "ERROR: wf must be larger than w0.\n";
        return 1;
    }

    // Safety: number of steps is non-positive
    if (nsteps <= 0)
    {
        std::cerr << "ERROR: 'nsteps' must be positive.\n";
        return 1;
    }

    // Output data file  
    std::ofstream output("data/moments_l" + std::to_string(l_plot) + ".dat");
    if (!output)
    {
        throw std::runtime_error("ERROR: Cannot open output file.");
    }
    output << std::scientific << std::setprecision(8); // setting number of digits in output file

    //==========================================//
    //===== EQUILIBRIUM INITIAL CONDITIONS =====//
    //==========================================//

    KiTable ki_table("data/Ki_table.dat");

    for (int i = nMin; i <= nMax; ++i)
    {
        for (int j = 0; j <= lMax; ++j)
        {
            rho[idx(i, j, lMax, nMin)] =
                (j == 0) ? equilibrium_moments(i, mass, initialTemp, alpha0, ki_table) : 0.0;
        }
    }

    double n0_initial = rho[idx(1, 0, lMax, nMin)];
    double e0_initial = rho[idx(2, 0, lMax, nMin)];
    double P0_initial = n0_initial * initialTemp; // classical ideal gas pressure
    double s0_initial = (e0_initial + P0_initial - alpha0 * initialTemp * n0_initial) / initialTemp;

    double tau_R0 = relaxation_time(eta_over_s, e0_initial, P0_initial, s0_initial);
    double t0 = w0 * tau_R0;
    double tau = std::log(t0);

    double dw = (wf - w0) / nsteps;
    double hmax = 1e-3;
    int max_substeps = 100000;

    //==========================================//
    //===== WRITE DATA IN OUTPUT DATA FILE =====//
    //==========================================//

    write_observables(output, tau, eta_over_s, mass, rho,
                      l_plot, nMin, nMax, lMax, ki_table);

    for (int step = 1; step <= nsteps; ++step)
    {
        double w_target = w0 + step * dw;

        int substeps = 0;

        while (scaled_time(tau, eta_over_s, mass, rho,
                           nMin, lMax, ki_table) < w_target)
        {
            if (substeps >= max_substeps)
            {
                throw std::runtime_error("ERROR: failed to reach target t/tau_R.");
            }

            double w_now = scaled_time(tau, eta_over_s, mass, rho,
                                       nMin, lMax, ki_table);

            double h_to_target = std::log(w_target / w_now);
            double h = std::min(hmax, h_to_target);

            if (h <= 1e-12)
            {
                throw std::runtime_error("ERROR: time step became too small.");
            }

            rk4_step(tau, eta_over_s, mass, rho, h,
                     nMin, nMax, lMax, ki_table);

            substeps++;
        }

        write_observables(output, tau, eta_over_s, mass, rho,
                          l_plot, nMin, nMax, lMax, ki_table);
    }

    return 0;
}


// version solving for t and plotting t/tau_R
// #include "rk4.hpp"
// #include "rhs.hpp"
// #include "KiTable.hpp"
// #include "utils.hpp"
// #include "constants.hpp"

// #include <cmath>
// #include <fstream>
// #include <iostream>
// #include <stdexcept>
// #include <string>
// #include <vector>
// #include <iomanip>

// using constants::GeVfm;

// void write_observables(std::ofstream& file,
//                        double tau,
//                        double eta_over_s,
//                        double mass,
//                        const std::vector<double>& rho,
//                        int l,
//                        int nmin,
//                        int nmax,
//                        int lmax,
//                        const KiTable& ki_table)
// {
//     double t = std::exp(tau);

//     double n0 = rho[idx(1, 0, lmax, nmin)];
//     double e0 = rho[idx(2, 0, lmax, nmin)];

//     double T = ki_table.T_Landau(mass, e0 / n0);
//     double alpha = ki_table.alpha_Landau(mass / T, T, n0);
//     double P0 = n0 * T; // classical ideal gas pressure
//     double s0 = (e0 + P0 - alpha * T * n0) / T;

//     double tau_R = 5.0 * eta_over_s * s0 / (e0 + P0);

//     file << t/tau_R;

//     for (int n = nmin; n <= nmax; ++n)
//     {
//         double rho_nl = rho[idx(n, l, lmax, nmin)];
//         double rho_eq_n0 = equilibrium_moments(n, mass, T, alpha, ki_table);

//         file << "\t" << std::abs(rho_nl) / rho_eq_n0;
//     }

//     file << "\n";
// }

// int main(int argc, char* argv[])
// {
//     if (argc != 5)
//     {
//         std::cerr << "Usage: " << argv[0] << " <l> <ti> <tf> <nsteps>\n";
//         return 1;
//     }

//     // Gather arguments parsed 
//     int l_plot = std::stoi(argv[1]);
//     double ti  = std::stod(argv[2]);
//     double tf  = std::stod(argv[3]);
//     int nsteps = std::stoi(argv[4]);

//     //==========================================//
//     //============= PHYSICAL INPUT =============//
//     //==========================================//

//     // Truncation
//     int nMin = 0;
//     int nMax = 5;
//     int lMax = 4;

//     int N = (nMax - nMin + 1) * (lMax + 1); // define dimension of the rho vector
//     std::vector<double> rho(N, 0.0); // allocate vector of corresponding size

//     // System properties
//     double eta_over_s = 0.5;
//     double mass = 0.135 * GeVfm;

//     // Initial conditions
//     double T0 = 0.2 * GeVfm;
//     double alpha0 = 0.0;
//     double tau0 = std::log(ti);
//     double tauf = std::log(tf);
//     double tau = tau0;
//     double h = (tauf - tau0) / nsteps;

//     //==========================================//
//     //================= GUARDS =================//
//     //==========================================//

//     // Check if chosen l is within the truncation
//     if (l_plot < 0 || l_plot > lMax)
//     {
//         std::cerr << "ERROR: l must satisfy 0 <= l <= " << lMax << "\n";
//         return 1;
//     }

//     // Check if initial and final times are positive
//     if (ti <= 0.0 || tf <= 0.0)
//     {
//         std::cerr << "ERROR: ti and tf must be positive.\n";
//         return 1;
//     }

//     // Safety: final time is smaller than initial time
//     if (tf <= ti)
//     {
//         std::cerr << "ERROR: tf must be larger than ti.\n";
//         return 1;
//     }

//     // Safety: number of steps is non-positive
//     if (nsteps <= 0)
//     {
//         std::cerr << "ERROR: 'nsteps' must be positive.\n";
//         return 1;
//     }

//     // Output data file  
//     std::ofstream output("data/moments_l" + std::to_string(l_plot) + ".dat");
//     if (!output)
//     {
//         throw std::runtime_error("ERROR: Cannot open output file.");
//     }
//     output << std::scientific << std::setprecision(8); // setting number of digits in output file

//     //==========================================//
//     //===== EQUILIBRIUM INITIAL CONDITIONS =====//
//     //==========================================//

//     KiTable ki_table("data/Ki_table.dat");

//     for (int i = nMin; i <= nMax; ++i)
//     {
//         for (int j = 0; j <= lMax; ++j)
//         {
//             rho[idx(i, j, lMax, nMin)] =
//                 (j == 0) ? equilibrium_moments(i, mass, T0, alpha0, ki_table) : 0.0;
//         }
//     }

//     //==========================================//
//     //===== WRITE DATA IN OUTPUT DATA FILE =====//
//     //==========================================//

//     for (int step = 0; step <= nsteps; ++step)
//     {
//         write_observables(output, tau, eta_over_s, mass, rho,
//                           l_plot, nMin, nMax, lMax, ki_table);

//         if (step < nsteps)
//         {
//             rk4_step(tau, eta_over_s, mass, rho, h,
//                      nMin, nMax, lMax, ki_table);
//         }
//     }

//     return 0;
// }

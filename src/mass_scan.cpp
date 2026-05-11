#include "rk4.hpp"
#include "rhs.hpp"
#include "KiTable.hpp"
#include "utils.hpp"
#include "constants.hpp"
#include "thermodynamics.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using constants::GeVfm;

int main(int argc, char* argv[])
{
    if (argc != 13)
    {

        std::cerr << "Usage: " << argv[0]
                  << " <tau0> <tauf> <nsteps> <mmin_GeV> <mmax_GeV> <nmass>"
                  << " <mode> <nFirst> <nLast> <lFirst> <lLast> <output_file>\n";
        return 1;
    }

    //==========================================//
    //============ PARSING ARGUMENTS ===========//
    //==========================================//

    double tau0      = std::stod(argv[1]);
    double tauf      = std::stod(argv[2]);
    int nsteps       = std::stoi(argv[3]);    
    double mmin_GeV  = std::stod(argv[4]);
    double mmax_GeV  = std::stod(argv[5]);
    int nmass        = std::stoi(argv[6]);
    std::string mode = argv[7];
    int nFirst       = std::stoi(argv[8]);
    int nLast        = std::stoi(argv[9]);
    int lFirst       = std::stoi(argv[10]);
    int lLast        = std::stoi(argv[11]);
    std::string output_filename = argv[12];

    //==========================================//
    //============ INPUT PARAMETERS ============//
    //==========================================//
    
    int nMin = 0;
    int nMax = 40;
    int lMax = 20;

    double eta_over_s = 1.0;
    double initialTemp = 10.0 * GeVfm;
    double alpha0 = 0.0;

    double tau_initial = std::log(tau0);
    double tau_final = std::log(tauf);
    double h = (tau_final - tau_initial) / nsteps;

    //==========================================//
    //================= GUARDS =================//
    //==========================================//

    if (tau0 <= 0.0 || tauf <= tau0 || nsteps <= 0)
        throw std::runtime_error("ERROR: invalid time interval or number of steps.");

    if (mmin_GeV < 0.0 || mmax_GeV < mmin_GeV || nmass <= 0)
        throw std::runtime_error("ERROR: invalid mass scan range.");

    if (mode != "n" && mode != "l")
        throw std::runtime_error("ERROR: mode must be either 'n' or 'l'.");

    if (nFirst < nMin || nLast > nMax || nLast < nFirst)
        throw std::runtime_error("ERROR: requested n range is outside the truncation.");

    if (lFirst < 0 || lLast > lMax || lLast < lFirst)
        throw std::runtime_error("ERROR: requested l range is outside the truncation.");

    if (mode == "n" && lFirst != lLast)
        throw std::runtime_error("ERROR: mode 'n' requires lFirst = lLast.");

    if (mode == "l" && nFirst != nLast)
        throw std::runtime_error("ERROR: mode 'l' requires nFirst = nLast.");


    //==========================================//
    //============= INSERT HEADERS =============//
    //==========================================//

    std::ofstream output(output_filename);
    if (!output)
        throw std::runtime_error("ERROR: Cannot open output file.");

    output << std::scientific << std::setprecision(8);

    // Mode
    output << "# mode " << mode << "\n";

    // Values of n
    output << "# n_values";
    if (mode == "n")
        for (int n = nFirst; n <= nLast; ++n) output << "\t" << n;
    else
        for (int l = lFirst; l <= lLast; ++l) output << "\t" << nFirst;
    output << "\n";

    // Values of l
    output << "# l_values";
    if (mode == "n")
        for (int n = nFirst; n <= nLast; ++n) output << "\t" << lFirst;
    else
        for (int l = lFirst; l <= lLast; ++l) output << "\t" << l;
    output << "\n";

    // Label columns of actual table
    output << "# columns\tmass_GeV";

    if (mode == "n")
        for (int n = nFirst; n <= nLast; ++n)
            output << "\t" << "abs_rho_" << n << "_" << lFirst << "_over_rhoeq_" << n << "_0";
    else
        for (int l = lFirst; l <= lLast; ++l)
            output << "\t" << "abs_rho_" << nFirst << "_" << l << "_over_rhoeq_" << nFirst << "_0";
    output << "\n";


    for (int imass = 0; imass < nmass; ++imass)
    {
        double mass_GeV = (nmass == 1)
            ? mmin_GeV
            : mmin_GeV + (mmax_GeV - mmin_GeV) * imass / (nmass - 1);

        double mass = mass_GeV * GeVfm;

        output << mass_GeV;
        
        //==========================================//
        //===== EQUILIBRIUM INITIAL CONDITIONS =====//
        //==========================================//
        
        int N = (nMax - nMin + 1) * (lMax + 1); // Compute length of the rho array from truncation
        std::vector<double> rho(N, 0.0);
        
        KiTable ki_table("data/Ki_table.dat"); // Create object KiTable

        for (int i = nMin; i <= nMax; ++i)
        {
            for (int j = 0; j <= lMax; ++j)
            {
                rho[idx(i, j, lMax, nMin)] =
                    (j == 0) ? equilibrium_moments(i, mass, initialTemp, alpha0, ki_table) : 0.0;
            }
        }

        double tau = tau_initial;

        for (int step = 0; step < nsteps; ++step)
        {
            rk4_step(tau, eta_over_s, mass, rho, h,
                     nMin, nMax, lMax, ki_table);
        }


        //==========================================//
        //==== EXPORT MOMENTS BASED ON MODE RUN ====//
        //==========================================//

        // Compute thermodynamic quantities
        double n0 = rho[idx(1, 0, lMax, nMin)];
        double e0 = rho[idx(2, 0, lMax, nMin)];
        thermodynamics thermo = compute_thermo(mass, eta_over_s, n0, e0, ki_table);

        if (mode == "n")
        {
            for (int n = nFirst; n <= nLast; ++n)
            {
                double rho_nl = rho[idx(n, lFirst, lMax, nMin)];
                double rho_eq_n0 = equilibrium_moments(n, mass, thermo.T, thermo.alpha, ki_table);
                output << "\t" << std::abs(rho_nl) / rho_eq_n0;
            }
        }
        else
        {
            for (int l = lFirst; l <= lLast; ++l)
            {
                double rho_nl = rho[idx(nFirst, l, lMax, nMin)];
                double rho_eq_n0 = equilibrium_moments(nFirst, mass, thermo.T, thermo.alpha, ki_table);
                output << "\t" << std::abs(rho_nl) / rho_eq_n0;
            }
        }

        output << "\n";
    }

    return 0;
}

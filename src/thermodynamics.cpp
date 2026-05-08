#include "thermodynamics.hpp"
#include "KiTable.hpp"
#include "utils.hpp"
#include "constants.hpp"

#include <cmath>
#include <stdexcept>

using constants::pi;

thermodynamics compute_thermo(double mass,
                              double eta_over_s,
                              double n0,
                              double e0,
                              const KiTable& ki_table)
{
    if (n0 <= 0.0)
        throw std::runtime_error("ERROR: particle density must be positive.");

    if (e0 <= 0.0)
        throw std::runtime_error("ERROR: energy density must be positive.");

    thermodynamics thermo;

    //==========================================//
    //======= LANDAU MATCHING CONDITIONS =======//
    //==========================================//

    // Massless gas
    if (mass < 1e-12)
    {
        thermo.T = e0 / (3.0 * n0);
        thermo.alpha = std::log((pi * pi * n0) / (std::pow(thermo.T, 3)));
    }
    // Massive gas
    else
    {
        thermo.T = ki_table.T_Landau(mass, e0 / n0);
        thermo.alpha = ki_table.alpha_Landau(mass / thermo.T, thermo.T, n0);
    }

    thermo.P0 = n0 * thermo.T; // classical ideal gas pressure
    thermo.s0 = (e0 + thermo.P0 - thermo.alpha * thermo.T * n0) / thermo.T;

    if (e0 + thermo.P0 <= 0.0)
        throw std::runtime_error("ERROR: e0 + P0 must be positive.");

    if (thermo.s0 <= 0.0)
        throw std::runtime_error("ERROR: entropy density became non-positive.");

    thermo.tau_R = 5.0 * eta_over_s * thermo.s0 / (e0 + thermo.P0);

    return thermo;
}

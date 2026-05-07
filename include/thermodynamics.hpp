#ifndef THERMODYNAMICS_HPP
#define THERMODYNAMICS_HPP

class KiTable;

struct thermodynamics
{
    double T;
    double alpha;
    double P0;
    double s0;
    double tau_R;
};

thermodynamics compute_thermo(double mass,
                              double eta_over_s,
                              double n0,
                              double e0,
                              const KiTable& ki_table);

#endif

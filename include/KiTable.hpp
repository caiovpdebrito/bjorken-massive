#ifndef KITABLE_HPP
#define KITABLE_HPP

#include <vector>
#include <string>
#include <cstddef>

class KiTable
{
private:

    std::vector<double> x_values;
    std::vector<int> n_values;
    std::vector<std::vector<double>> Ki;
    std::vector<double> k2k1_ratio; // stores the values of ĸ_2/ĸ_1 that will be used to compute T 
    size_t num_xpoints;   // number of x grid points

public:

    // Constructor
    KiTable(const std::string& filename);

    // Access functions
    const std::vector<double>& get_x() const { return x_values; }
    const std::vector<int>& get_n() const { return n_values; }
    const std::vector<std::vector<double>>& get_Ki() const { return Ki; }
    const std::vector<double>& get_k2k1_ratio() const { return k2k1_ratio; }
    int get_num_xpoints() const { return num_xpoints; }

    // Access single value
    double Ki_value(int n_index, int x_index) const { return Ki[n_index][x_index]; }

    // Compute ĸ_n(z) for values of z outside the list via interpolation
    double kappa(int n, double z) const;

    // Interpolate the ratio k2/k1
    double k2_over_k1(double z) const;

    // Compute temperature via root search
    double T_Landau(double mass, double e_over_n) const; 

    // Compute thermal potential
    double alpha_Landau(double z, double T, double n0) const;
};

#endif
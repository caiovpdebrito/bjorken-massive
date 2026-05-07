#include "KiTable.hpp"
#include "constants.hpp"
#include "utils.hpp"
#include <algorithm>   // std::find
#include <cmath>       // std::log
#include <fstream>     // std::ifstream
#include <sstream>     // std::stringstream
#include <stdexcept>   // std::runtime_error

using constants::pi;

// class KiTable
// {
// private:

//     std::vector<double> x_values;
//     std::vector<int> n_values;
//     std::vector<std::vector<double>> Ki;
//     std::vector<double> k2k1_ratio; // stores the values of ĸ_2/ĸ_1 that will be used to compute T 
//     size_t num_xpoints;   // number of x grid points

// public:

//     // Constructor
//     KiTable(const std::string& filename);

//     // Access functions
//     const std::vector<double>& get_x() const { return x_values; }
//     const std::vector<int>& get_n() const { return n_values; }
//     const std::vector<std::vector<double>>& get_Ki() const { return Ki; }
//     const std::vector<double>& get_k2k1_ratio() const { return k2k1_ratio; }
//     int get_num_xpoints() const { return num_xpoints; }

//     // Access single value
//     double Ki_value(int n_index, int x_index) const { return Ki[n_index][x_index]; }

//     // Compute ĸ_n(z) for values of z outside the list via interpolation
//     double kappa(int n, double z) const;

//     // Interpolate the ratio k2/k1
//     double k2_over_k1(double z) const;

//     // Compute temperature via root search
//     double T_Landau(double mass, double e_over_n) const; 

//     // Compute thermal potential
//     double alpha_Landau(double z, double T, double n0) const;
// };

//==========================================//
//=============== CONSTRUCTOR ==============//
//==========================================//

// Construct class object -- table containing the values of Bickley-Naylor functions
KiTable::KiTable(const std::string& filename)
{
    std::ifstream file(filename);

    if (!file)
        throw std::runtime_error("ERROR: Cannot open table file.");

    
    // Read header
    std::string line;
    std::getline(file, line);
    std::stringstream header(line);

    std::string token;
    header >> token; // skip "x"

    // Fill vector of values of n within the truncation
    int n;
    while (header >> n)
        n_values.push_back(n);

    int n_count = n_values.size();

    Ki.resize(n_count);

    // Read rows
    while (std::getline(file, line))
    {
        std::stringstream row(line);

        double x;
        row >> x;

        x_values.push_back(x);

        for (int i = 0; i < n_count; i++)
        {
            double value;
            row >> value;

            Ki[i].push_back(value);
        }
    }

    // Store number of x points
    num_xpoints = x_values.size();

    //==========================================//
    //============ COMPUTE ĸ_2/ĸ_1 =============//
    //==========================================//
    
    k2k1_ratio.resize(x_values.size());

    // Find iterators
    auto it1 = std::find(n_values.begin(), n_values.end(), 1);
    auto it2 = std::find(n_values.begin(), n_values.end(), 2);

    // Safety check
    if(it1 == n_values.end() || it2 == n_values.end())
        throw std::runtime_error("ERROR: Need n=1 and n=2 in table.");

    // Convert iterator -> index
    size_t n1 = std::distance(n_values.begin(), it1);
    size_t n2 = std::distance(n_values.begin(), it2);

    for (size_t i = 0; i < x_values.size(); i++)
    {
        k2k1_ratio[i] = Ki[n2][i] / Ki[n1][i];
    }

    file.close();

}

//==========================================//
//============== INTERPOLATION =============//
//==========================================//

// Compute the value of kappa_n(x) for a value of x not contained in the list via interpolation
double KiTable::kappa(int n, double x) const
{
    // Find index of requested n
    auto it = std::find(n_values.begin(), n_values.end(), n);

    // Check if requested n is in the list
    if (it == n_values.end()) // if the requested element is not found, the function 'find' returns the last element in the list
        throw std::runtime_error("ERROR: Requested n not in table.");

    size_t n_index = std::distance(n_values.begin(), it);

    // Check bounds
    if (x < x_values.front() || x > x_values.back())
        throw std::runtime_error("ERROR: x outside interpolation range.");

    // Find interval in x grid
    for (size_t i = 0; i < x_values.size() - 1; i++)
    {
        if (x >= x_values[i] && x <= x_values[i+1]) // check the index i for which x_{i} ≤ x ≤ x_{i+j}
        {
            double x0 = x_values[i];
            double x1 = x_values[i+1];

            double y0 = Ki[n_index][i];
            double y1 = Ki[n_index][i+1];

            double t = (x - x0) / (x1 - x0);

            return y0 + t * (y1 - y0);
        }
    }

    throw std::runtime_error("ERROR: Interpolation failed.");
}

// Interpolate the ratio ĸ2/ĸ1 
// [convenient since this function is constantly used to compute the temperature (Landau matching conditions)]
double KiTable::k2_over_k1(double z) const
{
    
    // for (size_t i = 0; i < x_values.size() - 1; i++)
    // {
    //     if (z >= x_values[i] && z <= x_values[i+1])
    //     {
    //         double z0 = x_values[i];
    //         double z1 = x_values[i+1];

    //         double r0 = k2k1_ratio[i];
    //         double r1 = k2k1_ratio[i+1];

    //         double t = (z - z0)/(z1 - z0);

    //         return r0 + t*(r1 - r0);
    //     }
    // }
    return kappa(2, z) / kappa(1, z);

    throw std::runtime_error("ERROR: value of z outside table.");
}

//==========================================//
//======= LANDAU MATCHING CONDITIONS =======//
//==========================================//

// Temperature
double KiTable::T_Landau(double mass, double e_over_n) const
{
    
    // Initial search values on the edges
    double xmin = x_values.front();
    double xmax = x_values.back();

    // Torelance value and maximum number of iterations
    double tol = 1e-6;
    int max_iter = 200;
    
    // Implementation of the bissection method
    int iter = 0;
    while((xmax - xmin) > tol && iter < max_iter)
    {
        double xmed = 0.5*(xmin + xmax); // compute midpoint of the current interval

        // Evaluate function at the and in the middle of the interval
        double fL = k2_over_k1(xmin) - xmin * e_over_n / mass;
        double fM = k2_over_k1(xmed) - xmed * e_over_n / mass;

        if(fL * fM <= 0.0) // check whether there is a zero between the two ends of the interval
            xmax = xmed; // if so, move the UPPER end of the interval to the midpoint and restart
        else
            xmin = xmed; // if not, move the LOWER end of the interval to the midpoint and restart

        iter++;
    }

    double z = 0.5*(xmin + xmax);

    return mass / z;
}

// Thermal potential
double KiTable::alpha_Landau(double z, double T, double n0) const
{
    return std::log((2*pi*pi) * (n0/intpower(T,3)) * 1/(kappa(1, z)));
}

//==========================================//
//============== TESTING AREA ==============//
//==========================================//

// int main()
// {
//     KiTable table("Ki_table.dat");

//     double value = table.k2_over_k1(2.46);

//     std::cout << "Interpolated kappa_n = " << value << std::endl;

//     // Testing zone
    
//     const auto& x = table.get_x();
//     const auto& n = table.get_n();
//     size_t numPoints = table.get_num_xpoints();

//     double temperature = table.T_Landau(0.15, 1);
//     double alpha = table.alpha_Landau(0.1, 0.15, 0.2);

//     // std::cout << "First x value: " << x[0] << std::endl;
//     // std::cout << "Number of points: " << numPoints << std::endl;
//     // std::cout << "First n value: " << n[0] << std::endl;
//     std::cout << "T = " << temperature << std::endl
//               << "a = " << alpha << std::endl;
    
//     return 0;
// }
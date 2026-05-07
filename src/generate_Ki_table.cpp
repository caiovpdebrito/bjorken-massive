#include "bickleynaylor.hpp"

#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char* argv[])
{
    if (argc != 6 && argc != 7)
    {
        std::cerr << "Usage: " << argv[0]
                  << " <nmin> <nmax> <xmin> <xmax> <npoints> [filename]\n";
        return 1;
    }

    // Gather arguments parsed
    int nMin = std::stoi(argv[1]);
    int nMax = std::stoi(argv[2]);
    double xMin = std::stod(argv[3]);
    double xMax = std::stod(argv[4]);
    int numPoints = std::stoi(argv[5]);

    std::string filename = "data/Ki_table.dat";

    if (argc == 7)
    {
        filename = argv[6];
    }

    //==========================================//
    //================= GUARDS =================//
    //==========================================//

    if (nMax < nMin)
    {
        std::cerr << "ERROR: nmax must be greater than or equal to nmin.\n";
        return 1;
    }

    if (xMin <= 0.0 || xMax <= 0.0)
    {
        std::cerr << "ERROR: xmin and xmax must be positive.\n";
        return 1;
    }

    if (xMax <= xMin)
    {
        std::cerr << "ERROR: xmax must be larger than xmin.\n";
        return 1;
    }

    if (numPoints <= 0)
    {
        std::cerr << "ERROR: npoints must be positive.\n";
        return 1;
    }

    write_Ki_table(filename,
                   nMin,
                   nMax,
                   xMin,
                   xMax,
                   numPoints);

    std::cout << "Table written to " << filename << "\n";

    return 0;
}

#include "rk4.hpp"
#include "rhs.hpp"
#include "KiTable.hpp"
#include "utils.hpp"
#include "constants.hpp"

#include <iostream>
#include <vector>
#include <cmath>

using constants::pi;
using constants::GeVfm;

int main(){

    // Truncation parameters
    int Imin = 0;
    int Imax = 5;
    int Jmax = 3;

    // System properties
    double eta_over_s = 1;
    double mass = 0.135*GeVfm; // gass of pions (mass in [fm^{-1}])

    // Initial conditions
    double initialTemp = 0.2*GeVfm; // initial temperature in fm^{-1}
    double alpha0 = 0; // initial thermal potential 
    double t0 = 1e-5; // initial time in fm
    double tau = std::log(t0);
    double h = 0.01;
    int nsteps = 1200;

    int N = (Imax - Imin + 1) * (Jmax + 1); // size of the vector based on the truncation
    std::vector<double> rho(N, 0.0); 

    //==========================================//
    //============ COMPUTE \kappa_n ============//
    //==========================================//

    KiTable ki_table("data/Ki_table.dat");

    // Equilibrium initial conditions
    for(int i = Imin; i <= Imax; ++i){
        for(int j = 0; j <= Jmax; ++j){
            rho[idx(i, j, Jmax, Imin)] = (j == 0) ? equilibrium_moments(i, mass, initialTemp, alpha0, ki_table) : 0; 
            // std::cout << rho[idx(i, j, Jmax, Imin)] << std::endl;
        }
    }

    for(int step=0; step<nsteps; ++step){
        rk4_step(tau, eta_over_s, mass, rho, h, Imin, Imax, Jmax, ki_table);

        if(step % 10 == 0){
            double t = std::exp(tau);
            //std::cout << "t=" << t << " rho[2,0]=" << rho[idx(2, 0, Jmax, Imin)] << std::endl;
            // std::cout << "t = " << t << " \t rho = " << rho[idx(1, 1, Jmax, Imin)] << std::endl;
            std::cout << rho[idx(1, 1, Jmax, Imin)] << std::endl;
            std::cout << rho[idx(2, 1, Jmax, Imin)] << std::endl;
            std::cout << rho[idx(3, 1, Jmax, Imin)] << std::endl;
            std::cout << rho[idx(1, 0, Jmax, Imin)] << std::endl;
            std::cout << rho[idx(2, 0, Jmax, Imin)] << std::endl;
        }
    }

    return 0;
}

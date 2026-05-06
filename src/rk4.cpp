#include "rk4.hpp"
#include "rhs.hpp"
#include "utils.hpp"
#include "KiTable.hpp"
#include <vector>
#include <cmath>
#include <iostream>

void rk4_step(double& tau,
              double eta_over_s,
              double mass,
              std::vector<double>& rho,
              double dt,
              int nmin, int nmax, int lmax,
              const KiTable& BickleyNaylor)
{
    std::vector<double> k1, k2, k3, k4, tmp; // allocate the RK4 coefficients and values of rho at each time step

    rhs(tau, eta_over_s, mass, rho, k1, nmin, nmax, lmax, BickleyNaylor);

    tmp = rho;
    for(size_t n=0; n<rho.size(); ++n) tmp[n] += 0.5*dt*k1[n];
    rhs(tau+0.5*dt, eta_over_s, mass, tmp, k2, nmin, nmax, lmax, BickleyNaylor);

    tmp = rho;
    for(size_t n=0; n<rho.size(); ++n) tmp[n] += 0.5*dt*k2[n];
    rhs(tau+0.5*dt, eta_over_s, mass, tmp, k3, nmin, nmax, lmax, BickleyNaylor);

    tmp = rho;
    for(size_t n=0; n<rho.size(); ++n) tmp[n] += dt*k3[n];
    rhs(tau+dt, eta_over_s, mass, tmp, k4, nmin, nmax, lmax, BickleyNaylor);

    for(size_t n=0; n<rho.size(); ++n)
        rho[n] += dt*(k1[n]+2*k2[n]+2*k3[n]+k4[n])/6.0;

    tau += dt;
}

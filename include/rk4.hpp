#ifndef RK4_HPP
#define RK4_HPP

#include <vector>

class KiTable;

void rk4_step(double& tau,
              double eta_over_s,
              double mass,
              std::vector<double>& rho,
              double dt,
              int nmin, int nmax, int lmax,
              const KiTable& BickleyNaylor);

#endif
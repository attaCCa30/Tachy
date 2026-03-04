#include <iostream>
#include <vector>
#include <cmath>

#include "parallel/MPIContext.hpp"
#include "bg/BackgroundState.hpp"
#include "bg/BackgroundStepper.hpp"
#include "pert/QEOM.hpp"
#include "driver/CoupledDriver.hpp"
#include "model/PaperPotential.hpp"

int main(int argc, char** argv) {
  MPIContext mpi(&argc, &argv);

  // ---- Parameters (edit here for now; later move to JSON config)
  PaperPotentialParams pp;
  pp.Mp = 1.0;
  pp.Vc = 1e-10;
  pp.alpha = 0.005;
  pp.A = 1.0;
  pp.B = 4.0;

  pp.m = 1e-6;
  pp.lambda = 1e-12;
  pp.g = 1e-3;
  pp.phi_SBP = 1.0;

  PaperPotential pot(pp);

  const double Mp = pp.Mp;
  BackgroundStepper bg_step(Mp, pot);
  QEOM q_eom(Mp, pot);

  // ---- k grid (log-spaced)
  const std::size_t Nk = 4096;
  const double kmin = 1e-5;
  const double kmax = 1e1;
  std::vector<double> k(Nk);
  for (std::size_t i=0; i<Nk; ++i) {
    const double f = double(i)/double(Nk-1);
    k[i] = kmin * std::pow(kmax/kmin, f);
  }
  const double dlnk = std::log(kmax/kmin)/double(Nk-1);

  // ---- initial background (choose consistent with your paper setup)
  BackgroundState bg0;
  bg0.N = 0.0;
  bg0.phi[0] = 10.0; // phi0
  bg0.phi[1] = 0.0;  // chi0
  bg0.phip[0] = 0.0; // dphi/dN
  bg0.phip[1] = 0.0; // dchi/dN
  bg0.eps = 0.0;
  bg0.H = std::sqrt(bg_step.H2(bg0, BackreactionSource{}));

  CoupledDriver driver(mpi, std::move(k), bg_step, q_eom, dlnk);

  driver.init_modes(bg0);

  // ---- time stepping
  const double dN = 1e-3;
  const int nSteps = 50000;
  const int outputEvery = 2000;

  if (mpi.rank==0) {
    std::cout << "Tachy: coupled background + Q^I modes with Hartree backreaction\n";
    std::cout << "MPI ranks=" << mpi.size << "\n";
  }

  driver.run(bg0, dN, nSteps, outputEvery, true);

  return 0;
}

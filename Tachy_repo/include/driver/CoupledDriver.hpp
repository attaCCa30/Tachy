#pragma once
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>

#include "parallel/MPIContext.hpp"
#include "parallel/KDecompose.hpp"

#include "bg/BackgroundState.hpp"
#include "bg/BackgroundStepper.hpp"
#include "backreaction/Backreaction.hpp"
#include "backreaction/VarianceBR.hpp"

#include "pert/ModeState.hpp"
#include "pert/InitialConditions.hpp"
#include "pert/QEOM.hpp"

#include "obs/Curvature.hpp"

struct CoupledDriver {
  MPIContext& mpi;
  KDecomp decomp;

  std::vector<double> k;
  std::vector<ModeMatrixState> modes_local;

  BackgroundStepper bg_step;
  QEOM q_eom;
  VarianceBR br_model;
  double dlnk;

  explicit CoupledDriver(MPIContext& mpi_,
                         std::vector<double> k_,
                         BackgroundStepper bg_step_,
                         QEOM q_eom_,
                         double dlnk_)
    : mpi(mpi_),
      decomp(KDecomp::block(k_.size(), mpi.rank, mpi.size)),
      k(std::move(k_)),
      modes_local(decomp.nloc()),
      bg_step(bg_step_), q_eom(q_eom_), dlnk(dlnk_) {}

  void init_modes(const BackgroundState& bg) {
    for (std::size_t j=0; j<modes_local.size(); ++j) {
      const std::size_t i = decomp.i0 + j;
      modes_local[j] = BD_IC(k[i], bg);
    }
  }

  void run(BackgroundState bg, double dN, int nSteps,
           int outputEvery,
           bool write_csv=true) {

    BackreactionSource br_global{};
    br_global.var_phi = 0.0;
    br_global.var_chi = 0.0;

    std::ofstream bg_out, spec_out;
    if (mpi.rank==0 && write_csv) {
      bg_out.open("output_background.csv");
      bg_out << "N,phi,chi,phip,chip,H,eps,var_phi,var_chi\n";
      spec_out.open("output_spectrum.csv");
      spec_out << "k,P_R\n";
    }

    // For a lightweight "final spectrum" output, we will gather P_R(k) at the end.
    // (Implementation: rank0 will gather P_R arrays via MPI_Gatherv in a later iteration.
    //  For now, we only compute and write spectrum on rank0 for the k-block that lives there.)

    for (int s=0; s<nSteps; ++s) {
      // Compute phipp at current step for C matrix
      Vec2 dphi_dN{}, dphip_dN{};
      bg_step.rhs(bg, br_global, dphi_dN, dphip_dN);
      const Vec2 phipp = dphip_dN;

      // Step modes (local block) with OpenMP
      #pragma omp parallel for schedule(static)
      for (std::size_t j=0; j<modes_local.size(); ++j) {
        const std::size_t i = decomp.i0 + j;
        modes_local[j] = q_eom.advance(modes_local[j], k[i], bg, br_global, phipp, dN);
      }

      // Local variances
      const BackreactionSource br_loc = br_model.compute_local(k, decomp.i0, modes_local, dlnk);

      // Global reduce
      br_global.var_phi = mpi.allreduce_sum(br_loc.var_phi);
      br_global.var_chi = mpi.allreduce_sum(br_loc.var_chi);

      // Step background with global backreaction
      bg = bg_step.advance(bg, br_global, dN);

      // periodic output on rank0
      if (mpi.rank==0 && outputEvery>0 && (s % outputEvery)==0) {
        std::cout << "step " << s
                  << " N=" << bg.N
                  << " eps=" << bg.eps
                  << " H=" << bg.H
                  << " var_chi=" << br_global.var_chi
                  << "\n";
        if (write_csv) {
          bg_out << bg.N << ","
                 << bg.phi[0] << "," << bg.phi[1] << ","
                 << bg.phip[0] << "," << bg.phip[1] << ","
                 << bg.H << "," << bg.eps << ","
                 << br_global.var_phi << "," << br_global.var_chi
                 << "\n";
        }
      }
    }

    // At the end, write partial spectrum (only for rank0's block) for quick plotting
    if (mpi.rank==0 && write_csv) {
      for (std::size_t j=0; j<modes_local.size(); ++j) {
        const std::size_t i = decomp.i0 + j;
        const double PR = P_R(k[i], bg, modes_local[j]);
        spec_out << k[i] << "," << PR << "\n";
      }
      bg_out.close();
      spec_out.close();
      std::cout << "Wrote output_background.csv and output_spectrum.csv (rank0 block only).\n";
      std::cout << "Next upgrade: MPI_Gatherv to assemble full spectrum on rank0.\n";
    }
  }
};

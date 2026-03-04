#pragma once
#include <vector>
#include <cmath>
#include "backreaction/Backreaction.hpp"
#include "pert/ModeState.hpp"

// Computes variances from the two basis solutions.
// For log-spaced k grid: integral d ln k.
struct VarianceBR {
  bool subtract_initial = false;
  // If subtract_initial, user should provide initial variances separately.
  double init_var_phi = 0.0;
  double init_var_chi = 0.0;

  BackreactionSource compute_local(const std::vector<double>& k,
                                  std::size_t i0,
                                  const std::vector<ModeMatrixState>& modes_local,
                                  double dlnk) const {
    BackreactionSource br{};
    const double norm = 1.0/(2.0*M_PI*M_PI);

    for (std::size_t j=0; j<modes_local.size(); ++j) {
      const std::size_t i = i0 + j;
      const double kk = k[i];
      const auto& MM = modes_local[j];

      // Sum over basis columns
      double Qphi2 = 0.0, Qchi2 = 0.0;
      for (int a=0; a<2; ++a) {
        const auto& m = MM.col[a];
        Qphi2 += m.Q.re[0]*m.Q.re[0] + m.Q.im[0]*m.Q.im[0];
        Qchi2 += m.Q.re[1]*m.Q.re[1] + m.Q.im[1]*m.Q.im[1];
      }

      br.var_phi += norm * dlnk * (kk*kk*kk) * Qphi2;
      br.var_chi += norm * dlnk * (kk*kk*kk) * Qchi2;
    }

    if (subtract_initial) {
      br.var_phi -= init_var_phi;
      br.var_chi -= init_var_chi;
      br.var_phi = std::max(0.0, br.var_phi);
      br.var_chi = std::max(0.0, br.var_chi);
    }

    return br;
  }
};

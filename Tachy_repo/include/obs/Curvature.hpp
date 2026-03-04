#pragma once
#include <cmath>
#include "bg/BackgroundState.hpp"
#include "pert/ModeState.hpp"

inline void R_components(const BackgroundState& bg, const ModeState& m,
                         double& Rre, double& Rim) {
  const double denom = bg.phip[0]*bg.phip[0] + bg.phip[1]*bg.phip[1];
  if (denom <= 0.0) { Rre=0.0; Rim=0.0; return; }

  const double num_re = bg.phip[0]*m.Q.re[0] + bg.phip[1]*m.Q.re[1];
  const double num_im = bg.phip[0]*m.Q.im[0] + bg.phip[1]*m.Q.im[1];

  Rre = num_re/denom;
  Rim = num_im/denom;
}

inline double R_abs2_from_matrix(const BackgroundState& bg, const ModeMatrixState& M) {
  // Quantum power: sum over basis columns |R_a|^2
  double sum = 0.0;
  for (int a=0; a<2; ++a) {
    double Rre=0.0, Rim=0.0;
    R_components(bg, M.col[a], Rre, Rim);
    sum += Rre*Rre + Rim*Rim;
  }
  return sum;
}

inline double P_R(double k, const BackgroundState& bg, const ModeMatrixState& M) {
  return (k*k*k)/(2.0*M_PI*M_PI) * R_abs2_from_matrix(bg, M);
}

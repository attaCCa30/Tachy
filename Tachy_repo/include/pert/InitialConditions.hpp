#pragma once
#include <cmath>
#include "pert/ModeState.hpp"
#include "bg/BackgroundState.hpp"

inline ModeState BD_IC_basis(double k, const BackgroundState& bg, int basis_field /*0 or 1*/) {
  ModeState m{};
  const double a = std::exp(bg.N);
  const double x = k/(a*bg.H); // k/(aH)
  const double amp = 1.0/std::sqrt(2.0*k);

  // basis vector e_I
  m.Q.re[basis_field] = amp;
  m.Q.im[basis_field] = 0.0;

  // Q_N = -i x Q  -> Re=0, Im=-x*Re
  m.P.re[basis_field] = 0.0;
  m.P.im[basis_field] = -x*amp;

  return m;
}

inline ModeMatrixState BD_IC(double k, const BackgroundState& bg) {
  ModeMatrixState M{};
  M.col[0] = BD_IC_basis(k, bg, 0);
  M.col[1] = BD_IC_basis(k, bg, 1);
  return M;
}

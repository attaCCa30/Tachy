#pragma once
#include <cmath>
#include "bg/BackgroundState.hpp"
#include "bg/Derived.hpp"
#include "model/Potential.hpp"
#include "backreaction/Backreaction.hpp"
#include "pert/ModeState.hpp"

struct QEOM {
  double Mp;
  const Potential& pot;
  explicit QEOM(double Mp_, const Potential& pot_) : Mp(Mp_), pot(pot_) {}

  Mat2 build_C_over_H2(const BackgroundState& bg,
                       const BackreactionSource& br,
                       const Vec2& phipp) const {
    const Mat2 V2 = pot.d2V_eff(bg.phi, br);
    const double H = bg.H;
    const double eps = bg.eps;

    // H'/H = -eps in N-time for canonical scalar field in FRW
    const double Hp_over_H = -eps;

    const Vec2& u = bg.phip;

    // X = H u u^T
    Mat2 X{};
    for (int i=0;i<NF;++i) for (int j=0;j<NF;++j) X[i][j] = H * u[i]*u[j];

    // X' = H' u u^T + H u' u^T + H u (u')^T
    Mat2 Xp{};
    for (int i=0;i<NF;++i) for (int j=0;j<NF;++j) {
      const double term1 = (Hp_over_H*H) * u[i]*u[j];
      const double term2 = H * phipp[i]*u[j];
      const double term3 = H * u[i]*phipp[j];
      Xp[i][j] = term1 + term2 + term3;
    }

    // D = (1/H) [ 3X + X' ]
    Mat2 D{};
    for (int i=0;i<NF;++i) for (int j=0;j<NF;++j) D[i][j] = (3.0*X[i][j] + Xp[i][j]) / H;

    // C/H^2 = V2/H^2 - D/H^2
    Mat2 C_over_H2{};
    const double H2 = H*H;
    for (int i=0;i<NF;++i) for (int j=0;j<NF;++j) {
      C_over_H2[i][j] = V2[i][j]/H2 - D[i][j]/H2;
    }
    return C_over_H2;
  }

  ModeState advance_one(const ModeState& m, double k, const BackgroundState& bg,
                        const BackreactionSource& br,
                        const Vec2& phipp, double dN) const {
    auto rhs = [&](const ModeState& s, ModeState& ds) {
      const double a = std::exp(bg.N);
      const double x = k/(a*bg.H);  // k/(aH)
      const double eps = bg.eps;
      const Mat2 M = build_C_over_H2(bg, br, phipp);

      ds.Q.re = s.P.re; ds.Q.im = s.P.im;

      for (int I=0; I<NF; ++I) {
        ds.P.re[I] = -(3.0 - eps)*s.P.re[I];
        ds.P.im[I] = -(3.0 - eps)*s.P.im[I];

        double sum_re = x*x * s.Q.re[I];
        double sum_im = x*x * s.Q.im[I];
        for (int J=0; J<NF; ++J) {
          sum_re += M[I][J]*s.Q.re[J];
          sum_im += M[I][J]*s.Q.im[J];
        }
        ds.P.re[I] += -sum_re;
        ds.P.im[I] += -sum_im;
      }
    };

    ModeState k1{}, k2{}, mid = m;
    rhs(m, k1);

    for (int I=0; I<NF; ++I) {
      mid.Q.re[I] += 0.5*dN*k1.Q.re[I];
      mid.Q.im[I] += 0.5*dN*k1.Q.im[I];
      mid.P.re[I] += 0.5*dN*k1.P.re[I];
      mid.P.im[I] += 0.5*dN*k1.P.im[I];
    }
    rhs(mid, k2);

    ModeState out = m;
    for (int I=0; I<NF; ++I) {
      out.Q.re[I] += dN*k2.Q.re[I];
      out.Q.im[I] += dN*k2.Q.im[I];
      out.P.re[I] += dN*k2.P.re[I];
      out.P.im[I] += dN*k2.P.im[I];
    }
    return out;
  }

  ModeMatrixState advance(const ModeMatrixState& M, double k, const BackgroundState& bg,
                          const BackreactionSource& br,
                          const Vec2& phipp, double dN) const {
    ModeMatrixState out = M;
    out.col[0] = advance_one(M.col[0], k, bg, br, phipp, dN);
    out.col[1] = advance_one(M.col[1], k, bg, br, phipp, dN);
    return out;
  }
};

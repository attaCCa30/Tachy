#pragma once
#include <cmath>
#include "bg/BackgroundState.hpp"
#include "bg/Derived.hpp"
#include "model/Potential.hpp"
#include "backreaction/Backreaction.hpp"

struct BackgroundStepper {
  double Mp;
  const Potential& pot;

  explicit BackgroundStepper(double Mp_, const Potential& pot_) : Mp(Mp_), pot(pot_) {}

  double H2(const BackgroundState& bg, const BackreactionSource& br) const {
    const double Veff = pot.V_eff(bg.phi, br);
    const double denom = 3.0*Mp*Mp - kineticN(bg);
    return Veff/denom;
  }

  void rhs(const BackgroundState& bg, const BackreactionSource& br,
           Vec2& dphi_dN, Vec2& dphip_dN) const {
    const double eps = epsilon(bg, Mp);
    const double h2  = H2(bg, br);
    const Vec2 grad  = pot.dV_eff(bg.phi, br);

    dphi_dN[0] = bg.phip[0];
    dphi_dN[1] = bg.phip[1];

    dphip_dN[0] = -(3.0 - eps)*bg.phip[0] - grad[0]/h2;
    dphip_dN[1] = -(3.0 - eps)*bg.phip[1] - grad[1]/h2;
  }

  BackgroundState advance(const BackgroundState& bg, const BackreactionSource& br, double dN) const {
    Vec2 k1_phi{}, k1_phip{};
    rhs(bg, br, k1_phi, k1_phip);

    BackgroundState mid = bg;
    mid.N += 0.5*dN;
    for (int I=0; I<NF; ++I) {
      mid.phi[I]  += 0.5*dN*k1_phi[I];
      mid.phip[I] += 0.5*dN*k1_phip[I];
    }

    Vec2 k2_phi{}, k2_phip{};
    rhs(mid, br, k2_phi, k2_phip);

    BackgroundState out = bg;
    out.N += dN;
    for (int I=0; I<NF; ++I) {
      out.phi[I]  += dN*k2_phi[I];
      out.phip[I] += dN*k2_phip[I];
    }

    out.eps = epsilon(out, Mp);
    out.H   = std::sqrt(H2(out, br));
    return out;
  }
};

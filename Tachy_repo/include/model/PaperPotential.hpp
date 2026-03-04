#pragma once
#include <cmath>
#include <algorithm>
#include "model/Potential.hpp"

/*
  Implements the paper potential "as-is" (two-field):
    V(phi) = Vc [1 + U(phi)]
    U(phi) = -1/2 * (phi^2/Mp^2) * ( B - A (1 + alpha ln(phi/Mp))^2 )

    V(chi) = -1/2 m^2 chi^2 + lambda chi^4

    V_int(phi,chi) = 1/2 g^2 chi^2 (phi - phi_SBP)^2

  Notes:
  - The paper uses m_Pl; in code we use Mp passed in parameters.
  - log(phi/Mp) requires phi>0. We use log(max(|phi|, tiny)/Mp) and keep derivatives with sign(phi).
    If you intend phi to remain positive, set phi0>0 and avoid crossing 0.
*/

struct PaperPotentialParams {
  double Mp = 1.0;   // reduced Planck mass in your units
  double Vc = 1e-10; // overall scale
  double alpha = 0.005;
  double A = 1.0;
  double B = 4.0;

  double m = 1e-6;       // hilltop mass parameter
  double lambda = 1e-12; // hilltop quartic
  double g = 1e-3;       // interaction coupling 
  double phi_SBP = 1.0;  // symmetry breaking point value
};

struct PaperPotential final : public Potential {
  PaperPotentialParams p;
  explicit PaperPotential(PaperPotentialParams pp) : p(pp) {}

  // Helpers
  inline double safe_abs(double x) const {
    return std::max(std::abs(x), 1e-300);
  }

  inline double ell(double phi) const {
    // log(|phi|/Mp)
    return std::log(safe_abs(phi)/p.Mp);
  }

  inline double one_plus_alpha_log(double phi) const {
    return 1.0 + p.alpha * ell(phi);
  }

  // U(phi)
  inline double U(double phi) const {
    const double x = phi / p.Mp;
    const double t = one_plus_alpha_log(phi);
    const double inside = p.B - p.A * (t*t);
    return -0.5 * (x*x) * inside;
  }

  double V(const Vec2& f) const override {
    const double phi = f[0], chi = f[1];
    const double Vphi = p.Vc * (1.0 + U(phi));
    const double Vchi = -0.5*p.m*p.m*chi*chi + p.lambda*chi*chi*chi*chi;
    const double dphi = (phi - p.phi_SBP);
    const double Vint = 0.5*p.g*p.g*chi*chi*dphi*dphi;
    return Vphi + Vchi + Vint;
  }

  Vec2 dV(const Vec2& f) const override {
    const double phi = f[0], chi = f[1];

    // dU/dphi
    // U = -1/2 (phi^2/Mp^2) (B - A t^2), t=1+α ln(|phi|/Mp)
    const double Mp = p.Mp;
    const double x = phi / Mp;
    const double t = one_plus_alpha_log(phi);
    const double inside = p.B - p.A*(t*t);

    // d/dphi (phi^2/Mp^2) = 2phi/Mp^2
    const double d_x2 = 2.0*phi/(Mp*Mp);

    // dt/dphi = α * d/dphi ln|phi| = α/phi (with sign handled by |phi| derivative -> 1/phi)
    // For phi != 0, d ln|phi| / dphi = 1/phi
    const double dt = p.alpha / ( (phi==0.0) ? ( (phi>=0)?1e-300:-1e-300 ) : phi );

    // dinside/dphi = -A*2 t dt
    const double dinside = -p.A * 2.0 * t * dt;

    // U = -1/2 * x^2 * inside
    // dU = -1/2*(d(x^2)*inside + x^2*dinside)
    const double dU = -0.5*( d_x2*inside + (x*x)*dinside );

    const double dVphi = p.Vc * dU;

    // dVchi
    const double dVchi = -p.m*p.m*chi + 4.0*p.lambda*chi*chi*chi;

    // Interaction gradients
    const double dphi = (phi - p.phi_SBP);
    const double g2 = p.g*p.g;
    const double dVint_dphi = g2*chi*chi*dphi;          // d/dphi 0.5 g2 chi^2 dphi^2 = g2 chi^2 dphi
    const double dVint_dchi = g2*chi*dphi*dphi;         // d/dchi ... = g2 chi dphi^2

    Vec2 g;
    g[0] = dVphi + dVint_dphi;
    g[1] = dVchi + dVint_dchi;
    return g;
  }

  Mat2 d2V(const Vec2& f) const override {
    const double phi = f[0], chi = f[1];
    const double Mp = p.Mp;
    const double x = phi / Mp;
    const double t = one_plus_alpha_log(phi);
    const double inside = p.B - p.A*(t*t);

    // dt/dphi = α/phi
    const double phi_safe = (phi==0.0) ? ((phi>=0)?1e-300:-1e-300) : phi;
    const double dt = p.alpha / phi_safe;
    // d2t/dphi2 = -α/phi^2
    const double d2t = -p.alpha / (phi_safe*phi_safe);

    // dinside = -2A t dt
    const double dinside = -2.0*p.A*t*dt;
    // d2inside = -2A[(dt)^2 + t*d2t]
    const double d2inside = -2.0*p.A*(dt*dt + t*d2t);

    // x^2 = phi^2/Mp^2
    const double x2 = x*x;
    const double dx2 = 2.0*phi/(Mp*Mp);
    const double d2x2 = 2.0/(Mp*Mp);

    // U = -1/2 x2 inside
    // dU = -1/2 (dx2 inside + x2 dinside)
    // d2U = -1/2 (d2x2 inside + 2 dx2 dinside + x2 d2inside)
    const double d2U = -0.5*( d2x2*inside + 2.0*dx2*dinside + x2*d2inside );

    const double Vphiph = p.Vc * d2U;

    // chi-chi from hilltop
    const double Vchich = -p.m*p.m + 12.0*p.lambda*chi*chi;

    // interaction hessian
    const double dphi = (phi - p.phi_SBP);
    const double g2 = p.g*p.g;
    const double Vint_phiphi = g2*chi*chi;              // d2/dphi2 0.5 g2 chi^2 dphi^2 = g2 chi^2
    const double Vint_chichi = g2*dphi*dphi;            // d2/dchi2 ... = g2 dphi^2
    const double Vint_phichi = 2.0*g2*chi*dphi;         // d2/dphi dchi ... = 2 g2 chi dphi

    Mat2 H{};
    H[0][0] = Vphiph + Vint_phiphi;
    H[1][1] = Vchich + Vint_chichi;
    H[0][1] = Vint_phichi;
    H[1][0] = Vint_phichi;
    return H;
  }

  // Hartree-effective:
  // - Use <chi^2> = chi^2 + var_chi, <(phi-phi*)^2> = (phi-phi*)^2 + var_phi in V_int
  // - For chi^4: <chi^4> = chi^4 + 6 chi^2 var_chi + 3 var_chi^2 (Gaussian)
  double V_eff(const Vec2& f, const BackreactionSource& br) const override {
    const double phi = f[0], chi = f[1];

    const double Vphi = p.Vc * (1.0 + U(phi));

    const double chi2eff = chi*chi + br.var_chi;
    const double Vchi = -0.5*p.m*p.m*chi2eff
                      + p.lambda*(chi*chi*chi*chi + 6.0*chi*chi*br.var_chi + 3.0*br.var_chi*br.var_chi);

    const double dphi = (phi - p.phi_SBP);
    const double phi2eff = dphi*dphi + br.var_phi;
    const double Vint = 0.5*p.g*p.g * chi2eff * phi2eff;

    return Vphi + Vchi + Vint;
  }

  Vec2 dV_eff(const Vec2& f, const BackreactionSource& br) const override {
    const double phi = f[0], chi = f[1];

    // keep V(phi) bare in mean-field (no Taylor var_phi term here)
    const Vec2 grad_bare = dV(f);

    const double dphi = (phi - p.phi_SBP);
    const double g2 = p.g*p.g;

    const double chi2eff = chi*chi + br.var_chi;
    const double phi2eff = dphi*dphi + br.var_phi;

    // Replace interaction gradients with Hartree-effective factors:
    const double dVint_dphi = g2*chi2eff*dphi;
    const double dVint_dchi = g2*phi2eff*chi;

    // Hilltop gradient with Hartree:
    // Vchi_eff = -1/2 m^2(chi^2+var) + lambda(chi^4 + 6 chi^2 var + 3 var^2)
    // d/dchi: -m^2 chi + 4λ chi^3 + 12 λ chi var
    const double dVchi_eff = -p.m*p.m*chi + 4.0*p.lambda*chi*chi*chi + 12.0*p.lambda*chi*br.var_chi;

    Vec2 g;
    // For phi: take bare phi part from grad_bare but replace interaction term
    g[0] = (grad_bare[0] - (g2*chi*chi*dphi)) + dVint_dphi;
    // For chi: replace entire chi gradient by Hartree version + interaction
    g[1] = dVchi_eff + dVint_dchi;
    return g;
  }

  Mat2 d2V_eff(const Vec2& f, const BackreactionSource& br) const override {
    const double phi = f[0], chi = f[1];
    const double dphi = (phi - p.phi_SBP);
    const double g2 = p.g*p.g;

    // Start from bare Hessian then patch terms:
    Mat2 H = d2V(f);

    const double chi2eff = chi*chi + br.var_chi;
    const double phi2eff = dphi*dphi + br.var_phi;

    // Interaction Hessian (Hartree):
    // d2/dphi2: g2 * chi2eff
    // d2/dchi2: g2 * phi2eff
    // d2/dphidchi: 2 g2 chi dphi (same form, with mean chi)
    const double Hphiphi_eff = g2 * chi2eff;
    const double Hchichi_eff = g2 * phi2eff;
    const double Hphichi_eff = 2.0*g2*chi*dphi;

    // Replace interaction pieces in H:
    // Bare interaction pieces were: g2*chi^2, g2*dphi^2, 2 g2 chi dphi
    // So update by subtracting bare and adding effective.
    H[0][0] += (Hphiphi_eff - g2*chi*chi);
    H[1][1] += (Hchichi_eff - g2*dphi*dphi);
    H[0][1] = Hphichi_eff;
    H[1][0] = Hphichi_eff;

    // Hilltop chi-chi with Hartree:
    // d2/dchi2: -m^2 + 12 λ chi^2 + 12 λ var_chi
    H[1][1] += 12.0*p.lambda*br.var_chi;

    return H;
  }
};

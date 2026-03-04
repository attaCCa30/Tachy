#pragma once
#include "core/Types.hpp"
#include "backreaction/Backreaction.hpp"

struct Potential {
  virtual ~Potential() = default;

  // Bare (mean-field) potential and derivatives
  virtual double V(const Vec2& f) const = 0;
  virtual Vec2 dV(const Vec2& f) const = 0;
  virtual Mat2 d2V(const Vec2& f) const = 0;

  // Hartree-effective versions. Default: no backreaction.
  virtual double V_eff(const Vec2& f, const BackreactionSource& br) const {
    (void)br; return V(f);
  }
  virtual Vec2 dV_eff(const Vec2& f, const BackreactionSource& br) const {
    (void)br; return dV(f);
  }
  virtual Mat2 d2V_eff(const Vec2& f, const BackreactionSource& br) const {
    (void)br; return d2V(f);
  }
};

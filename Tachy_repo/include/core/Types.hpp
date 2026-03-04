#pragma once
#include <array>

constexpr int NF = 2; // (phi, chi)

struct Vec2 {
  double v[NF]{0.0, 0.0};
  double& operator[](int i) { return v[i]; }
  double  operator[](int i) const { return v[i]; }
};

struct Mat2 {
  double m[NF][NF]{{0.0,0.0},{0.0,0.0}};
  double* operator[](int i) { return m[i]; }
  const double* operator[](int i) const { return m[i]; }
};

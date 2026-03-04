#pragma once
#include <cmath>
#include "bg/BackgroundState.hpp"

inline double kineticN(const BackgroundState& bg) {
  return 0.5*(bg.phip[0]*bg.phip[0] + bg.phip[1]*bg.phip[1]);
}
inline double epsilon(const BackgroundState& bg, double Mp) {
  return kineticN(bg)/(Mp*Mp);
}

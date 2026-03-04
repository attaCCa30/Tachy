#pragma once
#include "core/Types.hpp"

struct BackgroundState {
  double N = 0.0;
  Vec2 phi;
  Vec2 phip; // dphi/dN
  double H = 0.0;
  double eps = 0.0;
};

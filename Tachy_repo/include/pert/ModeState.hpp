#pragma once
#include "core/Types.hpp"

struct CplxVec2 {
  Vec2 re;
  Vec2 im;
};

struct ModeState {
  CplxVec2 Q; // Q^I
  CplxVec2 P; // dQ^I/dN
};

struct ModeMatrixState {
  // Two independent basis solutions (columns)
  ModeState col[2];
};

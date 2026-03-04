#pragma once
#include <cstddef>

struct KDecomp {
  std::size_t i0=0, i1=0; // [i0,i1)
  static KDecomp block(std::size_t Nk, int rank, int size) {
    KDecomp d;
    d.i0 = (static_cast<std::size_t>(rank) * Nk) / static_cast<std::size_t>(size);
    d.i1 = (static_cast<std::size_t>(rank+1) * Nk) / static_cast<std::size_t>(size);
    return d;
  }
  std::size_t nloc() const { return i1 - i0; }
};

# Tachy (linear perturbations + backreaction; Q^I core)

This repository is a starter C++ codebase to study the trapped-inflation mechanism with:
- 2-field homogeneous background evolution
- multi-field Mukhanov–Sasaki perturbations `Q^I_k` (flat-gauge / gauge-invariant core)
- Hartree-style backreaction from mode variances feeding the background
- hybrid MPI + OpenMP parallelism over `k` modes

## Build (MPI + OpenMP)
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DENABLE_MPI=ON -DENABLE_OPENMP=ON
cmake --build build -j
```

## Run (example: 8 MPI ranks x 16 threads = 128 cores)
```bash
export OMP_NUM_THREADS=16
export OMP_PROC_BIND=close
export OMP_PLACES=cores
mpirun -np 8 ./build/tachy
```

Output:
- prints periodic diagnostics to stdout
- writes `output_background.csv` and `output_spectrum.csv` on rank 0 (lightweight; plot in Python)

## Notes
- The paper potential is implemented in `include/model/PaperPotential.hpp`.
- Backreaction currently uses a simple "vacuum-subtracted" option hook but defaults to no subtraction.
  For production runs you should implement adiabatic subtraction or at least subtract the initial BD variance.

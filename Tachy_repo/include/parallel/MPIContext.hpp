#pragma once
#ifdef ENABLE_MPI
  #include <mpi.h>
#endif

struct MPIContext {
  int rank=0, size=1;

  MPIContext(int* argc, char*** argv) {
#ifdef ENABLE_MPI
    MPI_Init(argc, argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#else
    (void)argc; (void)argv;
#endif
  }

  ~MPIContext() {
#ifdef ENABLE_MPI
    MPI_Finalize();
#endif
  }

  double allreduce_sum(double local) const {
#ifdef ENABLE_MPI
    double global=0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global;
#else
    return local;
#endif
  }
};

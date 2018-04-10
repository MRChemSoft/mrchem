option(ENABLE_BLAS "Enable use of BLAS? (Requires a Fortran compiler)" OFF)
if(ENABLE_BLAS)
  include(FindBLAS)
  if(BLAS_FOUND)
    # Enable BLAS in Eigen
    add_definitions(-DEIGEN_USE_BLAS)
  endif(BLAS_FOUND)
endif(ENABLE_BLAS)

option(ENABLE_BLAS "Enable use of BLAS? (Requires a Fortran compiler)" OFF)
if(ENABLE_BLAS)
  include(FindBLAS)
  if(BLAS_FOUND)
    set(HAVE_BLAS 1)
    # Get symbol names
    enable_language(Fortran)
    include(FortranCInterface)
    FortranCInterface_HEADER(FCMangle.h
      MACRO_NAMESPACE "FC_"
      SYMBOL_NAMESPACE "FC_"
      SYMBOLS dgemm ddot)
  endif(BLAS_FOUND)
endif(ENABLE_BLAS)

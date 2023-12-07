I have built and ran three different sets of calculations

the ones denoted with only mpi are run with mrchem built from the main branch at MRChemsoft.

The ones with vecmult have used the vecmult branch.
specifically:

for `mpi_vecmult`:
Built and installed MRCPP from vecmult branch with omp and mpi
Used this MRCPP to build MRChem from vecmult branch with omp and mpi

for `omp_vecmult`:
built and installed MRCPP from vecmult with only omp
Used this MRCPP to build MRCHem from vecmult branch with only omp


The main branch mpi has issues with clearing the functions at the end of the 
microiteration loops. This has been addressed by Peter.

The vecmult branch has no issues with this in MPI and it doesnt seem to affect the omp procedure either. 
The vecmult branch is not updating the update for the microiterations.

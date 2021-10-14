# Change log

## Version 1.0.2 2021-10-13

### Changed

- Updated MRCPP to v1.4.0
- Updated Eigen to v3.4.0

### Fixed

- Compile error with Eigen-3.4 due to conversion std::array<double> -> Eigen::VectorXd

## Version 1.0.1 2020-12-04

### Changed

- Moved CI from Travis to GitHub Actions
- Abort on failure when reading MW guess

### Fixed

- Minor issues with the printed output
- Error in NMR tensor output (dia term printed transposed)
- Precision issue with CHK guess (now full world_prec is used)
- MPI dead-lock when mixing shared-mem potentials with orbital bank
- SAD guess for large systems (Eigen threading error in matrix multiply)
- CORE/SAD guess for heavy atoms (force refinement for sharp s-functions)

## Version 1.0.0 2020-10-28

### Added

- Installation instructions to README

### Changed

- Updated MRCPP to v1.3.6

### Fixed

- Faulty MPI_INTEGER type

## Version 1.0.0-alpha3 2020-10-10

### Added

- Dagger terms in exchange hessian
- Screening of exchange contributions

### Changed

- Updated MRCPP to v1.3.5
- Updated XCFun to v2.1.0
- Updated parselglossy to v0.7.0
- Removed runtest dependency
- Removed parselglossy as runtime dependency
- Encapsulated MRChem and MRCPP parallelizations

### Fixed

- Installation path for SAD files
- Misc CMake fixes to enable packaging

## Version 1.0.0-alpha2 2020-06-29

### Added

- New JSON output

### Changed

- Improved mrchem launcher script
- Improved documentation
- Input keywords for v1 defined and fixed
- Updated MRCPP to v1.2.0
- Updated XCFun to v2.0.1
- Updated nlohmann/json to v3.6.1

## Version 1.0.0-alpha1 2020-05-05

### Features

- Hartree-Fock
- Kohn-Sham DFT (LDA/GGA/hybrid)
- Restricted (closed-shell) and unrestricted
- External electric field
- Ground state energy
- Dipole moment
- Quadrupole moment
- Polarizability (static/dynamic)
- Magnetizability
- NMR shielding
- Density/orbital plots

## Version 0.2.2 2019-11-20

### Fixed

- Updated MRCPP to v1.0.2
- OpenMPI error with complex data types

## Version 0.2.1 2019-08-12

### Fixed

- Updated MRCPP to v1.0.1
- Eigen installation


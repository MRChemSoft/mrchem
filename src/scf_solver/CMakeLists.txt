target_sources(mrchem PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/Accelerator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/EnergyOptimizer.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/HelmholtzVector.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/KAIN.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/OrbitalOptimizer.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GroundStateSolver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/LinearResponseSolver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SCF.cpp
    )

target_include_directories(mrchem
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  )

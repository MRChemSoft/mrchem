target_sources(mrchem PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/AOBasis.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/AOContraction.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/AOFileReader.h
    ${CMAKE_CURRENT_SOURCE_DIR}/DaltonAOReader.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GTOFileReader.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Intgrl.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/OrbitalExp.cpp
  )

# Readers for HDF5-based GTO files (DIRAC, Serenity, VeloxChem) go here, e.g.
# if(ENABLE_HDF5)
#   target_sources(mrchem PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/VeloxChemAOReader.cpp)
# endif()

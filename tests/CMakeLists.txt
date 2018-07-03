# Prepare MRChem "library" for other executables
#set(MRChem_INCLUDE_DIR ${CMAKE_BINARY_DIR}/stage${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR})
#set(MRChem_LIBRARY_DIR ${CMAKE_BINARY_DIR}/stage${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBRARYDIR})
#add_library(MRChem INTERFACE)
#target_include_directories(MRChem INTERFACE ${MRChem_INCLUDE_DIR})
#target_link_libraries(MRChem INTERFACE ${MRChem_LIBRARY_DIR})

add_executable(mrchem-tests unit_tests.cpp)

add_subdirectory(qmfunctions)
add_subdirectory(qmoperators)

target_link_libraries(mrchem-tests
  PRIVATE
    mrchem
    Catch
    Eigen3::Eigen
  )

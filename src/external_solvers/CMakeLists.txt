target_sources(mrchem PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/ExternalSolver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ChemTensorSolver.cpp
)
target_link_libraries(mrchem PRIVATE chemtensor_lib)
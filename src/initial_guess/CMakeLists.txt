# Install sources in lib
target_sources(mrchem PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/core.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gto.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sad.cpp
)

#####################################
# Install initial guess executables #
#####################################

# Install core guess executable in bin
add_executable(core-guess core_guess.cpp)
target_link_libraries(core-guess
  PRIVATE
    mrchem
  )
install(
  TARGETS core-guess
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

# Install SAD guess executable in bin
add_executable(sad-guess sad_guess.cpp)
target_link_libraries(sad-guess
  PRIVATE
    mrchem
  )
install(TARGETS sad-guess
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Install GTO guess executable in bin
add_executable(gto-guess gto_guess.cpp)
target_link_libraries(gto-guess
  PRIVATE
    mrchem
  )
install(
  TARGETS
    gto-guess
  RUNTIME DESTINATION
    ${CMAKE_INSTALL_BINDIR}
  )

# Install driver for LSDalton calculations on stallo
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/lsdalton/lsdalton-guess.in lsdalton-guess)
install(
  PROGRAMS
    ${CMAKE_CURRENT_BINARY_DIR}/lsdalton-guess
  DESTINATION
    ${CMAKE_INSTALL_BINDIR}
  )

install(
  FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/lsdalton/lsdalton.py
  DESTINATION
    ${PYTHON_SITE_INSTALL_DIR}
  )

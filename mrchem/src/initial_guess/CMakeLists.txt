target_sources(mrchem
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/core.cpp
    ${CMAKE_CURRENT_LIST_DIR}/sad.cpp
    ${CMAKE_CURRENT_LIST_DIR}/gto.cpp
  )
# TODO These should be exported in the the configuration file

# <<< Core guess executable >>>
add_executable(core-guess core_guess.cpp)
target_link_libraries(core-guess
  PRIVATE
    mrchem
  )

install(
  TARGETS core-guess
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

# <<< SAD guess executable >>>
add_executable(sad-guess sad_guess.cpp)
target_link_libraries(sad-guess
  PRIVATE
    mrchem
  )

install(TARGETS sad-guess
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# <<< GTO guess executable >>>
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

# <<< Driver for LSDalton calculations on stallo >>>
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

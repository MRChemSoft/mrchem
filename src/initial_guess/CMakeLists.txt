target_sources(mrchem
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/chk.cpp
    ${CMAKE_CURRENT_LIST_DIR}/core.cpp
    ${CMAKE_CURRENT_LIST_DIR}/gto.cpp
    ${CMAKE_CURRENT_LIST_DIR}/mw.cpp
    ${CMAKE_CURRENT_LIST_DIR}/sad.cpp
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

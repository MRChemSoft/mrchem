# Try to find LibXC. We only enable the backend if we both (a) opt in with
# MRCHEM_ENABLE_LIBXC=ON and (b) actually find a working Libxc::xc target.

find_package(Libxc QUIET)

# Some distros provide an imported target Libxc::xc, some only variables.
if(TARGET Libxc::xc)
  set(HAVE_LIBXC TRUE)
  get_property(_loc TARGET Libxc::xc PROPERTY LOCATION)
  message(STATUS "Found LibXC: ${_loc}")
elseif(LIBXC_FOUND)
  # Compat path: construct an IMPORTED target from legacy variables.
  add_library(Libxc::xc UNKNOWN IMPORTED)
  set_target_properties(Libxc::xc PROPERTIES
    IMPORTED_LOCATION "${LIBXC_LIBRARIES}"
    INTERFACE_INCLUDE_DIRECTORIES "${LIBXC_INCLUDE_DIRS}"
  )
  set(HAVE_LIBXC TRUE)
  message(STATUS "Found LibXC (legacy variables).")
else()
  set(HAVE_LIBXC FALSE)
  message(STATUS "LibXC not found – LibXC backend will be disabled unless you install it.")
endif()

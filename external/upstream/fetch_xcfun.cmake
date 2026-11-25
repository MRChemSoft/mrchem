option(XCFUN_FIND_BEHAVIOUR "Find XCFun behavior" "default")
cpm_set_find_behaviour(${XCFUN_FIND_BEHAVIOUR})
set(CMAKE_BUILD_TYPE Release)
if(NOT XCFUN_OLD_PBE)
  add_compile_definitions(XCFUN_REF_PBEX_MU)
endif()
CPMAddPackage(
  NAME XCFun
  VERSION 2.1
  GIT_TAG v2.1.1
  GITHUB_REPOSITORY dftlibs/xcfun
  FIND_PACKAGE_ARGUMENTS "CONFIG NO_CMAKE_PATH NO_CMAKE_PACKAGE_REGISTRY NO_CMAKE_SYSTEM_PACKAGE_REGISTRY"
  OPTIONS
  "ENABLE_TESTALL OFF"
  "XCFUN_MAX_ORDER 3"
  "XCFUN_PYTHON_INTERFACE OFF"
  )
if(TARGET xcfun)
  target_compile_options(xcfun PRIVATE -w)
  if(XCFUN_OLD_PBE)
    message(STATUS "Compiling XCFun with old PBE parameters (different from LibXC)")
  else()
    message(STATUS "Compiling XCFun with new PBE parameters (same as LibXC)")
  endif()
else()
  message(STATUS "Using local xcfun: mrchem does not have control over which PBE parameters are used (might be different from LibXC)")
endif()
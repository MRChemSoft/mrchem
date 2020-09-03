find_package(nlohmann_json 3.6 CONFIG QUIET
  NO_CMAKE_PATH
  NO_CMAKE_SYSTEM_PACKAGE_REGISTRY
  )
if(TARGET nlohmann_json::nlohmann_json)
  get_target_property(
    _loc
    nlohmann_json::nlohmann_json
    INTERFACE_INCLUDE_DIRECTORIES
  )
  message(STATUS "Found nlohmann_json: ${_loc} (found version ${nlohmann_json_VERSION})")
else()
  message(STATUS "Suitable nlohmann_json could not be located: downloading and building nlohmann_json instead.")
  include(FetchContent)
  FetchContent_Declare(nlohmann_json_sources
    QUIET
    URL
      https://github.com/nlohmann/json/archive/v3.6.1.tar.gz
    )

  FetchContent_GetProperties(nlohmann_json_sources)

  set(JSON_BuildTests OFF CACHE BOOL "" FORCE)
  set(JSON_ImplicitConversions OFF CACHE BOOL "" FORCE)

  if(NOT nlohmann_json_sources_POPULATED)
    FetchContent_Populate(nlohmann_json_sources)

    add_subdirectory(
      ${nlohmann_json_sources_SOURCE_DIR}
      ${nlohmann_json_sources_BINARY_DIR}
      )
  endif()
endif()

if(NLOHMANN_JSON_FIND_BEHAVIOUR STREQUAL "default" OR NLOHMANN_JSON_FIND_BEHAVIOUR STREQUAL "onlylocal")
  find_package(nlohmann_json 3.12 CONFIG QUIET
    NO_CMAKE_PATH
    NO_CMAKE_PACKAGE_REGISTRY
    NO_CMAKE_SYSTEM_PACKAGE_REGISTRY
    )
endif()
if(TARGET nlohmann_json::nlohmann_json)
  get_target_property(_loc nlohmann_json::nlohmann_json INTERFACE_INCLUDE_DIRECTORIES)
  message(STATUS "Found nlohmann_json: ${_loc} (found version ${nlohmann_json_VERSION})")
elseif(NLOHMANN_JSON_FIND_BEHAVIOUR STREQUAL "default" OR NLOHMANN_JSON_FIND_BEHAVIOUR STREQUAL "onlyfetch")
  message(STATUS "Suitable nlohmann_json could not be located. Fetching and building!")
  include(FetchContent)
  FetchContent_Declare(nlohmann_json_sources
    QUIET
    URL
      https://github.com/nlohmann/json/archive/v3.12.0.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )

  set(JSON_BuildTests OFF CACHE BOOL "" FORCE)
  set(JSON_ImplicitConversions ON CACHE BOOL "" FORCE)

  FetchContent_MakeAvailable(nlohmann_json_sources)
else()
  message(FATAL_ERROR "No suitable nlohmann_json found or fetched. Aborting setup!")
endif()

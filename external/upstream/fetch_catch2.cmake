find_package(Catch2 2.7 CONFIG QUIET
  NO_CMAKE_PATH
  NO_CMAKE_PACKAGE_REGISTRY
  )

if(TARGET Catch2::Catch2)
  message(STATUS "Using Catch2: ${Catch2_SOURCE_DIR} (version ${Catch2_VERSION})")
else()
  message(STATUS "Suitable Catch2 could not be located. Fetching and building!")
  include(FetchContent)

  include(FetchContent)
  FetchContent_Declare(Catch2
    QUIET
    URL https://github.com/catchorg/Catch2/archive/v2.13.4.tar.gz
  )

  FetchContent_MakeAvailable(Catch2)
endif()

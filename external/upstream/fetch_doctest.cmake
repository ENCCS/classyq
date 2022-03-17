find_package(doctest 2.4 CONFIG QUIET
  NO_CMAKE_PATH
  NO_CMAKE_PACKAGE_REGISTRY
  )

if(TARGET doctest::doctest)
  message(STATUS "Using doctest: ${doctest_SOURCE_DIR} (version ${doctest_VERSION})")
else()
  message(STATUS "Suitable doctest could not be located. Fetching and building!")
  include(FetchContent)

  include(FetchContent)
  FetchContent_Declare(doctest
    QUIET
    URL https://github.com/doctest/doctest/archive/v2.4.8.tar.gz
  )

  FetchContent_MakeAvailable(doctest)
endif()

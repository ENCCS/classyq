find_package(autodiff 0.6 CONFIG QUIET
  NO_CMAKE_PATH
  NO_CMAKE_PACKAGE_REGISTRY
  )

if(TARGET autodiff::autodiff)
  message(STATUS "Found autodiff version ${autodiff_VERSION}")
else()
  message(STATUS "Suitable autodiff could not be located. Fetching!")
  include(FetchContent)

  FetchContent_Declare(autodiff
    QUIET
    URL
      https://github.com/autodiff/autodiff/archive/v0.6.6.zip
    )

  set(AUTODIFF_BUILD_TESTS    OFF CACHE BOOL "" FORCE)
  set(AUTODIFF_BUILD_PYTHON   OFF CACHE BOOL "" FORCE)
  set(AUTODIFF_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
  set(AUTODIFF_BUILD_DOCS     OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(autodiff)
endif()

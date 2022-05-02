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
    GIT_REPOSITORY
      https://github.com/autodiff/autodiff
    GIT_TAG
     58e1ca5afd20c5b3b85529729a5dac8fc87c2561  # main branch on 2022-05-02 
    #URL
    #  https://github.com/autodiff/autodiff/archive/v0.6.7.zip
    )

  set(AUTODIFF_BUILD_TESTS    OFF CACHE BOOL "" FORCE)
  set(AUTODIFF_BUILD_PYTHON   OFF CACHE BOOL "" FORCE)
  set(AUTODIFF_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
  set(AUTODIFF_BUILD_DOCS     OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(autodiff)
endif()

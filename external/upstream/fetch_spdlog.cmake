find_package(spdlog 1.9.2 CONFIG QUIET
  NO_CMAKE_PATH
  NO_CMAKE_PACKAGE_REGISTRY
  )

if(TARGET spdlog::spdlog)
  message(STATUS "Found spdlog version ${spdlog_VERSION}")
else()
  message(STATUS "Suitable spdlog could not be located. Fetching!")
  include(FetchContent)

  FetchContent_Declare(spdlog
    QUIET
    URL
      https://github.com/gabime/spdlog/archive/v1.9.2.zip
    )

  FetchContent_MakeAvailable(spdlog)
endif()

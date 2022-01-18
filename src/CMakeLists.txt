# we compile everything as object files
add_library(classyq-objs OBJECT LeopardiPartition.cpp Sphere.cpp)

set_target_properties(classyq-objs
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
  )

target_compile_features(classyq-objs
  PUBLIC
    cxx_std_17
  )

target_compile_definitions(classyq-objs
  INTERFACE
    $<INSTALL_INTERFACE:USING_ClassyQ>
    SPDLOG_FMT_EXTERNAL
    H5_USE_EIGEN
  )

if(WIN32)
  target_compile_definitions(classyq-objs
    PUBLIC
      _USE_MATH_DEFINES
      _CRT_SECURE_NO_WARNINGS
      _SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING
    )
endif()

target_include_directories(classyq-objs
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )

target_include_directories(classyq-objs
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}
  INTERFACE
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/ClassyQ>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/ClassyQ>
  )

target_link_libraries(classyq-objs
  PUBLIC
    Eigen3::Eigen
    fmt::fmt
    HighFive
    spdlog::spdlog
  )

#add_subdirectory(...)

add_executable(foo foo.cpp)

target_compile_features(foo
  PUBLIC
    cxx_std_17
  )

target_link_libraries(foo
  PRIVATE
    classyq-objs
  )

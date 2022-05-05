# we compile everything as object files
add_library(classyq-objs OBJECT)

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
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/ClassyQ>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/ClassyQ>
  )

foreach(_dir IN ITEMS cavity utils)
  add_subdirectory(${_dir})
endforeach()

target_link_libraries(classyq-objs
  PUBLIC
    Eigen3::Eigen
    fmt::fmt
    HighFive
    spdlog::spdlog
  )

foreach(exe IN ITEMS foo)
add_executable(${exe} ${exe}.cpp)

target_link_libraries(${exe}
  PRIVATE
    classyq-objs
  )
endforeach()

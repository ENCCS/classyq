# download Catch2 if not available on the system
include(${PROJECT_SOURCE_DIR}/external/upstream/fetch_catch2.cmake)

add_executable(unit_tests
  main.cpp
  add.cpp
  )

target_link_libraries(unit_tests
  PRIVATE
    classyq-objs
    Catch2::Catch2
  )

# automatic discovery of unit tests
list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/contrib)
include(Catch)
catch_discover_tests(unit_tests
  PROPERTIES
    LABELS "unit"
  )

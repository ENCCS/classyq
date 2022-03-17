# download doctest if not available on the system
include(${PROJECT_SOURCE_DIR}/external/upstream/fetch_doctest.cmake)

add_executable(unit_tests
  main.cpp
  )

target_link_libraries(unit_tests
  PRIVATE
    classyq-objs
    doctest::doctest
  )

foreach(_dir IN ITEMS cavity)
  add_subdirectory(${_dir})
endforeach()

# automatic discovery of unit tests
include(doctest)
doctest_discover_tests(unit_tests
  PROPERTIES
    LABELS "unit"
  )

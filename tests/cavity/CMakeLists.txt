target_sources(unit_tests
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/leopardi.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/tsless.cpp
  )

configure_file(leopardi.h5 leopardi.h5 COPYONLY)

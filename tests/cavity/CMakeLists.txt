target_sources(unit_tests
  PRIVATE
    leopardi.cpp
  )

configure_file(leopardi.h5 leopardi.h5 COPYONLY)

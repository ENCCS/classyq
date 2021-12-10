# Export compile commands for each file to JSON
# This is useful for static analysis tools and linters
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_EXTENSIONS OFF)

add_subdirectory(src)

enable_testing()
include(CTest)

# This must come last!!
add_subdirectory(tests)

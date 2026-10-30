cmake_minimum_required(VERSION 3.10)
project(reportable_shortcuts CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")
# Hardware popcount and, where available, BMI2 (for select within a word).
option(NATIVE "compile for this machine's processor (-march=native)" ON)
if(NATIVE)
  add_compile_options(-march=native)
endif()

# libsais (vendored, Apache 2.0), for suffix arrays of small inputs.
add_library(sais64 STATIC third_party/libsais/src/libsais64.c third_party/libsais/src/libsais.c)
target_include_directories(sais64 PUBLIC third_party/libsais/include)

add_library(shortcuts STATIC
  src/bitvector.cpp
  src/int_vector.cpp
  src/wavelet_tree.cpp
  src/text.cpp
  src/rlbwt.cpp
  src/shortcuts.cpp
  src/build.cpp
  src/search.cpp
  src/mems.cpp
)
target_include_directories(shortcuts PUBLIC include)
target_link_libraries(shortcuts PUBLIC sais64)

add_executable(build_index tools/build_index.cpp)
target_link_libraries(build_index PRIVATE shortcuts)

add_executable(query_index tools/query_index.cpp)
target_link_libraries(query_index PRIVATE shortcuts)

add_executable(test_bitvector tests/test_bitvector.cpp)
target_link_libraries(test_bitvector PRIVATE shortcuts)

add_executable(test_shortcuts tests/test_shortcuts.cpp)
target_link_libraries(test_shortcuts PRIVATE shortcuts)

add_executable(test_mems tests/test_mems.cpp)
target_link_libraries(test_mems PRIVATE shortcuts)

add_executable(query_mems tools/query_mems.cpp)
target_link_libraries(query_mems PRIVATE shortcuts)

enable_testing()
add_test(NAME test_bitvector COMMAND test_bitvector)
add_test(NAME test_shortcuts COMMAND test_shortcuts)
add_test(NAME test_mems COMMAND test_mems)

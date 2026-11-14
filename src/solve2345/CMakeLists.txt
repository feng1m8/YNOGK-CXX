cmake_minimum_required(VERSION 3.26)

project(
    solve2345
    LANGUAGES C
)

file(GENERATE OUTPUT ".gitignore" CONTENT "*")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_EXTENSIONS OFF)

file(
    REAL_PATH "../../extra/solve2345"
    SOLVE2345_CURRENT_SOURCE_DIR
)

add_library(
    solve2345 STATIC
    ${SOLVE2345_CURRENT_SOURCE_DIR}/src/quadratic.c
    ${SOLVE2345_CURRENT_SOURCE_DIR}/src/cubic.c
    ${SOLVE2345_CURRENT_SOURCE_DIR}/src/quartic.c
    ${SOLVE2345_CURRENT_SOURCE_DIR}/src/quintic.c
)

target_compile_definitions(solve2345 PRIVATE FP_PREC=53)

target_include_directories(
    solve2345 INTERFACE
    ${SOLVE2345_CURRENT_SOURCE_DIR}
)

add_library(
    solve2345l OBJECT
    ${SOLVE2345_CURRENT_SOURCE_DIR}/src/quadratic.c
    ${SOLVE2345_CURRENT_SOURCE_DIR}/src/cubic.c
    ${SOLVE2345_CURRENT_SOURCE_DIR}/src/quartic.c
    ${SOLVE2345_CURRENT_SOURCE_DIR}/src/quintic.c
)

target_compile_definitions(solve2345l PRIVATE FP_PREC=64)

target_compile_options(
    solve2345l PUBLIC
    -include preinclude.h
)

target_include_directories(
    solve2345l PUBLIC
    "."
    ${SOLVE2345_CURRENT_SOURCE_DIR}/src
)

target_link_libraries(solve2345 PRIVATE solve2345l)

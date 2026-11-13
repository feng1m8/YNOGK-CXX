add_library(
    solve2345 OBJECT
    src/solve2345/src/quadratic.c
    src/solve2345/src/cubic.c
    src/solve2345/src/quartic.c
    src/solve2345/src/quintic.c
)

target_compile_definitions(
    solve2345 PRIVATE
    FP_PREC=53
)

target_compile_options(
    solve2345 PRIVATE
    -include src/patch/solve2345.h
)

target_include_directories(
    solve2345 
    PRIVATE "."
    PRIVATE "src/solve2345/src"
    INTERFACE "src"
)

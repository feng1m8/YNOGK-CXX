cmake_minimum_required(VERSION 3.26)

set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project(
    "YNOGK-CXX"
    LANGUAGES CXX C
    VERSION "2023.1.5"
)

file(GENERATE OUTPUT ".gitignore" CONTENT "*")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_EXTENSIONS OFF)

add_library(
    "ynogk_cxx" OBJECT
    "YNOGK/common_src/Carlsons.c"
    "YNOGK/common_src/ellCarlsons.c"
    "YNOGK/common_src/ellfunctions.c"
    "YNOGK/common_src/particle.c"
    "YNOGK/common_src/rmu_tps.c"
    "YNOGK/common_src/root34.c"
    "src/radius.cxx"
    "src/theta.cxx"
    "src/impl.cxx"
    "src/particle.cxx"
)

target_include_directories(
    "ynogk_cxx" 
    PUBLIC "include" 
    PUBLIC "YNOGK/common_src"
)

find_package(
    Python
    COMPONENTS Interpreter Development.SABIModule NumPy
)

if(Python_FOUND)
    list(APPEND CMAKE_MODULE_PATH "${Python_SITEARCH}/cython_cmake/cmake")

    find_package(Cython MODULE)

    if(Cython_FOUND)
        include(UseCython)
        
        cython_transpile(
            "src/ynogk.pyx"
            LANGUAGE CXX
            OUTPUT_VARIABLE ynogk_pyx
        )

        Python_add_library(
            "ynogk" MODULE
            USE_SABI 3.8 WITH_SOABI
            ${ynogk_pyx}
        )

        target_link_libraries(
            "ynogk"
            PRIVATE ynogk_cxx 
            PRIVATE Python::NumPy
        )

        install(
            TARGETS "ynogk"
            DESTINATION "."
        )
    endif()
endif()

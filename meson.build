project(
    'YNOGK-CXX', 'cpp', 'c',
    version: '2023.1.5'
    meson_version: '>=1.3.0',
    license: 'GPL-3.0-or-later',
    license_files: [
        'LICENSE',
        'YNOGK/LICENSE',
    ],
    default_options: {
        'b_lto': true,
        'b_ndebug': 'if-release',
        'buildtype': 'release',
        'cpp_std': 'c++17',
        'c_std': 'c11',
    },
)

ynogk_cxx = static_library(
    'ynogk_cxx',
    'YNOGK/common_src/Carlsons.c',
    'YNOGK/common_src/ellCarlsons.c',
    'YNOGK/common_src/ellfunctions.c',
    'YNOGK/common_src/particle.c',
    'YNOGK/common_src/rmu_tps.c',
    'YNOGK/common_src/root34.c',
    'src/radius.cxx',
    'src/theta.cxx',
    'src/impl.cxx',
    'src/particle.cxx',
    include_directories: [
        'include',
        'YNOGK/common_src',
    ],
    build_by_default: false,
)

ynogk_cxx = declare_dependency(
    objects: [
        ynogk_cxx.extract_all_objects(),
    ],
    include_directories: [
        'include',
        'YNOGK/common_src',
    ],
)

python = import('python').find_installation(
    pure: false,
    modules: [
        'numpy',
        'cython',
    ],
    required: false,
)

if python.found()
    add_languages('cython')

    numpy = dependency('numpy')

    python.extension_module(
        'ynogk',
        'src/ynogk.pyx',
        dependencies: [
            ynogk_cxx,
            numpy,
        ],
        override_options: [
            'cython_language=cpp',
        ],
        limited_api: '3.8',
        install: true,
    )
endif

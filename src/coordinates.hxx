#ifndef YNOGK_CXX_COORDINATES_HXX
#define YNOGK_CXX_COORDINATES_HXX

#if defined(_MSC_VER) and !defined(_CRT_USE_C_COMPLEX_H)
#define _CRT_USE_C_COMPLEX_H
#ifndef _C_COMPLEX_T
#define _C_COMPLEX_T
typedef double _Complex _C_double_complex;
typedef float _Complex _C_float_complex;
typedef long double _Complex _C_ldouble_complex;
#endif
#endif

#include <complex.h>

extern "C"
{
#include "BLcoordinates_new.h"
#include "ellCarlsons.h"
}

namespace ynogk
{
    class Coordinates
    {
    public:
#define static
#include "BLcoordinates_new.c"
#undef static

    private:
        void exit(int) {}
    };
}

#endif

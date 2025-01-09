#ifndef YNOGK_CXX_COORDINATES_HXX
#define YNOGK_CXX_COORDINATES_HXX

#include <complex.h>

extern "C"
{
#include "ellCarlsons.h"
#include "BLcoordinates_new.h"
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

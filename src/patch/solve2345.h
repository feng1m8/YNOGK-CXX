#ifndef SOLVE2345_PATCH_H
#define SOLVE2345_PATCH_H

#include <stdlib.h>
#include <math.h>

#include "solve_utils.h"

#if FP_PREC == 53
#undef SCALB
#define SCALB scalbn

#undef EPS
#define EPS (1.0 / (1ull << 52))
#endif

#endif

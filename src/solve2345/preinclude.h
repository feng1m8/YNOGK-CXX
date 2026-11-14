#ifndef SOLVE2345_PREINCLUDE_H
#define SOLVE2345_PREINCLUDE_H

#include <float.h>
#include <math.h>
#include <stdlib.h>

#include "solve_utils.h"

#if FP_PREC == 53
#undef EPS
#define EPS (1.0 / (1ull << 52))
#undef SCALB
#define SCALB scalbn
#endif

#if FP_PREC == 64
#undef EPS
#define EPS (1.0 / (1ull << 63))
#undef SCALB
#define SCALB scalbnl

#if LDBL_MAX_EXP <= 4000
#undef EPS
#define EPS (1.0 / (1ull << 52))
#undef HYPOT
#define HYPOT hypotl
#undef FP_PREC
#define FP_PREC 53
#endif
#endif

#endif

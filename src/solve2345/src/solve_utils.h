#ifndef _SOLVE_UTILS_H
#define _SOLVE_UTILS_H

#define boolean int
#define false 0
#define true  1

#if FP_PREC == 53
#define EPS (1.0 / (1ul << 52))
#define FP_TYPE double
#define SQRT    sqrt
#define CBRT    cbrt
#define FABS    fabs
#define LOG     log
#define POW     pow
#define SCALB   scalb
#define FREXP   frexp
#define HYPOT   hypot
#define MIN_NORMAL DBL_MIN
#define MAX_VALUE  DBL_MAX
#endif

#if FP_PREC == 64
#define EPS (1.0 / (1ul << 63))
#define FP_TYPE long double
#define SQRT    sqrtl
#define CBRT    cbrtl
#define FABS    fabsl
#define LOG     logl
#define POW     powl
#define SCALB   scalbl
#define FREXP   frexpl
#define HYPOT   hypotl_impl
#define MIN_NORMAL LDBL_MIN
#define MAX_VALUE  LDBL_MAX
#endif

#if FP_PREC == 64
#define quadraticSolve quadraticSolvel
#define quadraticSolveMonic quadraticSolveMonicl
#define quadraticCoefs quadraticCoefsl

#define cubicSolve cubicSolvel
#define cubicSolveMonic cubicSolveMonicl
#define cubicSolveDepressed cubicSolveDepressedl
#define cubicSolveDepressedDominant cubicSolveDepressedDominantl
#define cubicCoefs cubicCoefsl

#define quarticSolve quarticSolvel
#define quarticSolveMonic quarticSolveMonicl
#define quarticSolveDepressed quarticSolveDepressedl
#define quarticSolveFullyDepressed quarticSolveFullyDepressedl
#define quarticSolveBiquadratic quarticSolveBiquadraticl
#define quarticCoefs quarticCoefsl

#define quinticSolve quinticSolvel
#define quinticSolveMonic quinticSolveMonicl
#define quinticSolveBringJerrardDominant quinticSolveBringJerrardDominantl
#define quinticSolveBringJerrard quinticSolveBringJerrardl
#define quinticCoefs quinticCoefsl
#endif

#ifndef FP_TYPE
#error FP_PREC must be set to either 53 or 64
#endif


#ifdef max
#undef max
#endif
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
      _a > _b ? _a : _b; })

#ifdef min
#undef min
#endif
#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
      _a < _b ? _a : _b; })

#ifdef abs
#undef abs
#endif
#define abs(a) \
   ({ __typeof__ (a) _a = (a); \
      _a < 0 ? -_a : _a; })



void quadraticSolve(FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE *zeror, FP_TYPE *zeroi);
void quadraticSolveMonic(FP_TYPE a, FP_TYPE b, FP_TYPE *zeror, FP_TYPE *zeroi);
void quadraticCoefs(FP_TYPE factor, int nReal, FP_TYPE a, FP_TYPE b, FP_TYPE *coefs);

void cubicSolve(FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE d, FP_TYPE *zeror, FP_TYPE *zeroi);
void cubicSolveMonic(FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE *zeror, FP_TYPE *zeroi);
void cubicSolveDepressed(FP_TYPE g, FP_TYPE h, FP_TYPE *zeror, FP_TYPE *zeroi);
FP_TYPE cubicSolveDepressedDominant(FP_TYPE g, FP_TYPE h);
void cubicCoefs(FP_TYPE factor, int nReal, FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE *coefs);

void quarticSolve(FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE d, FP_TYPE e, FP_TYPE *zeror, FP_TYPE *zeroi);
void quarticSolveMonic(FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE d, FP_TYPE *zeror, FP_TYPE *zeroi);
void quarticSolveDepressed(FP_TYPE b, FP_TYPE c, FP_TYPE d, FP_TYPE *zeror, FP_TYPE *zeroi);
void quarticSolveFullyDepressed(FP_TYPE c, FP_TYPE d, FP_TYPE *zeror, FP_TYPE *zeroi);
void quarticSolveBiquadratic(FP_TYPE a, FP_TYPE c, FP_TYPE e, FP_TYPE *zeror, FP_TYPE *zeroi);
void quarticCoefs(FP_TYPE factor, int nReal, FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE d, FP_TYPE *coefs);

void quinticSolve(FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE d, FP_TYPE e, FP_TYPE f,
                  FP_TYPE *zeror, FP_TYPE *zeroi);
void quinticSolveMonic(FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE d, FP_TYPE e,
                       FP_TYPE *zeror, FP_TYPE *zeroi);
FP_TYPE quinticSolveBringJerrardDominant(FP_TYPE g, FP_TYPE h);
void quinticSolveBringJerrard(FP_TYPE g, FP_TYPE h, FP_TYPE *zeror, FP_TYPE *zeroi);
void quinticCoefs(FP_TYPE factor, int nReal, FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE d, FP_TYPE e, FP_TYPE *coefs);

#endif   /* _SOLVE_UTILS_H */

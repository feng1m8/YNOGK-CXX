#ifndef _TEST53_H
#define _TEST53_H

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



static long double relErrorReal(int n, long double a, long double *zeror, long double *zeroi) {
    if (a == 0) {
        return 0;
    }

    long double err = LDBL_MAX;
    long double a1 = 1.0/abs(a);
    for (int i=0; i<n; i++) {
        long double e = a1*hypotl(zeror[i] - a, zeroi[i]);
        err = min(err, e);
    }
    return err;
}



static long double relErrorComplex(int n, long double a, long double b, long double *zeror, long double *zeroi) {
    if (a == 0 && b == 0) {
        return 0;
    }

    long double err = LDBL_MAX;
    long double a1 = 1.0/hypotl(a, b);
    for (int i=0; i<n; i++) {
        long double e = a1*hypotl(zeror[i] - a, zeroi[i] - b);
        err = min(err, e);
    }
    return err;
}


#endif
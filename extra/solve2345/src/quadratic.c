#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "solve_utils.h"

static inline FP_TYPE diff_of_products(FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE d);
static inline int getExponent(FP_TYPE a);


/**
 *
 * @author woelen
 */

/**
 * Solves a quadratic equation with real coefficients of the form a*x^2 + b*x + c == 0.
 * The complete solution is returned in the supplied arrays, in ascending order of the
 * real part and for a complex conjugate pair, the root with positive imaginary part
 * is given first.
 * @param a The coefficient of the quadratic term
 * @param b The coefficient of the linear term
 * @param c The constant term
 * @param zeror This is used for output. The first two elements of this array are
 * overwritten with the real part of the two zeros of the equation.
 * @param zeroi This is used for output. The first two elements of this array are
 * overwritten with the imaginary part of the two zeros of the equation.
 */    
void quadraticSolve(FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE *zeror, FP_TYPE *zeroi) {

    // If there is no constant term, then reduce the
    // equation to a linear equation. Assure that the
    // roots are ordered in ascending order.
    if (c == 0) {
        FP_TYPE root = -b/a;
        if (root < 0) {
            zeror[0] = root;
            zeror[1] = 0.0;
        }
        else {
            zeror[1] = root;
            zeror[0] = 0.0;
        }
        zeroi[0] = zeroi[1] = 0.0;
        return;
    }

    // Use b', such that we get the the equation a*x^2 + 2*b'*x + c == 0.
    // The formulae in terms of b' are somewhat simpler than the standard
    // ABC-formula for b.
    //    2*b' = b, so b' = b/2;
    FP_TYPE ba = b/2;

    // We use a scaled version of the discriminant to
    // avoid overflow or underflow. Overflow or underflow
    // may still occur, but only when the roots are close
    // to the representable range of a FP_TYPE (i.e. when
    // the magnitude is well over 10^300 or when the 
    // magnitude is well below 10^(-300)). If we want
    // complete protection against overflow or underflow,
    // then we would need a much more elaborate scaling
    // mechanism, but the cost in terms of complexity and
    // slowing down the code is too high to justify the
    // use of complete protection.
    // We scale by the maximum of the approximate binary
    // exponents of |b| and sqrt(|a*c|). By scaling by
    // a power of 2 we do not lose any precision.
    int exp = -max(getExponent(ba), (getExponent(a) + getExponent(c))/2);
    FP_TYPE D = diff_of_products(SCALB(ba, exp), SCALB(ba, exp), SCALB(a, exp), SCALB(c, exp));

    if (D < 0) {
        // Two complex conjugate roots. First the one with positive
        // imaginary part is stored and then the one with negative
        // imaginary part.
        FP_TYPE s = SCALB(SQRT(-D), -exp);
        zeror[0] = zeror[1] = -ba/a;
        zeroi[0] = ((a < 0) ? (-s) : s)/a;
        zeroi[1] = -zeroi[0];
    }
    else {
        // Two real roots. Store them in ascending order.
        FP_TYPE s = SCALB(SQRT(D), -exp);
        FP_TYPE bs = (ba < 0) ? (-ba + s) : (-ba - s);
        FP_TYPE root0 = bs/a;
        FP_TYPE root1 = c/bs;
        if (root0 < root1) {
            zeror[0] = root0;
            zeror[1] = root1;
        }
        else {
            zeror[1] = root0;
            zeror[0] = root1;
        }
        zeroi[0] = zeroi[1] = 0.0;
    }
}



/**
 * Solves a quadratic equation with real coefficients of the form x^2 + a*x + b == 0.
 * The complete solution is returned in the supplied arrays, in no particular order.
 * @param a The coefficient of the linear term
 * @param b The constant term
 * @param zeror This is used for output. The first two elements of this array are
 * overwritten with the real part of the two zeros of the equation.
 * @param zeroi This is used for output. The first two elements of this array are
 * overwritten with the imaginary part of the two zeros of the equation.
 */    
void quadraticSolveMonic(FP_TYPE a, FP_TYPE b, FP_TYPE *zeror, FP_TYPE *zeroi) {

    // If there is no constant term, then reduce the
    // equation to a linear equation. Assure that the
    // roots are ordered in ascending order.
    if (b == 0) {
        FP_TYPE root = -a;
        if (root < 0) {
            zeror[0] = root;
            zeror[1] = 0.0;
        }
        else {
            zeror[0] = 0.0;
            zeror[1] = root;
        }
        zeroi[0] = zeroi[1] = 0.0;
        return;
    }

    // Use a', such that we get the the equation x^2 + 2*a'*x + b == 0.
    //    2*a' = a, so a' = a/2;
    FP_TYPE aa = a/2;

    // We use a scaled version of the discriminant to
    // avoid overflow or underflow. Overflow or underflow
    // may still occur, but only when the roots are close
    // to the representable range of a FP_TYPE (i.e. when
    // the magnitude is well over 10^300 or when the 
    // magnitude is well below 10^(-300)). If we want
    // complete protection against overflow or underflow,
    // then we would need a much more elaborate scaling
    // mechanism, but the cost in terms of complexity and
    // slowing down the code is too high to justify the
    // use of complete protection.
    // We scale by the maximum of the approximate binary
    // exponents of |b| and sqrt(|a*c|). By scaling by
    // a power of 2 we do not lose any precision.
    int exp = -max(getExponent(aa), getExponent(b)/2);
#if FP_PREC == 53
    FP_TYPE D = fma(SCALB(aa, exp), SCALB(aa, exp), -SCALB(b, 2*exp));
#endif
#if FP_PREC == 64
    FP_TYPE D = SCALB(aa, exp) * SCALB(aa, exp) - SCALB(b, 2*exp);
#endif

    if (D < 0) {
        // Two complex conjugate roots. First the one with positive
        // imaginary part is stored and then the one with negative
        // imaginary part.
        FP_TYPE s = SCALB(SQRT(-D), -exp);
        zeror[0] = zeror[1] = -aa;
        zeroi[0] = s;
        zeroi[1] = -zeroi[0];
    }
    else {
        // Two real roots. Store them in ascending order.
        FP_TYPE s = SCALB(SQRT(D), -exp);
        FP_TYPE bs = (aa < 0) ? (-aa + s) : (-aa - s);
        FP_TYPE root0 = bs;
        FP_TYPE root1 = b/bs;
        if (root0 < root1) {
            zeror[0] = root0;
            zeror[1] = root1;
        }
        else {
            zeror[1] = root0;
            zeror[0] = root1;
        }
        zeroi[0] = zeroi[1] = 0.0;
    }
}


    
    
    
    
    
/**
 * Given a set of roots and the coefficient of the leading term, compute
 * all coefficients of the quadratic equation.
 * @param factor The coefficient of the leading term (the quadratic term).
 * @param nReal Tells whether the equation has two real roots, or zero real
 * roots. The next two parameters have meaning, depending on this parameter.
 * @param a If nReal equals 2, then this is interpreted as one of the real
 * roots of the equation. If nReal equals 0, then this is interpreted as the
 * real part of both complex conjugate roots.
 * @param b If nReal equals 2, then this is interpreted as one of the real
 * roots of the equation. If nReal equals 0, then this is interpreted as the
 * imaginary part of one of the complex conjugate roots, the other complex
 * root has imaginary part equal to -b.
 * @param coefs An output array, which is filled with the coefficients
 * of the quadratic. This array must have a length at least equal to 3. 
 * At index 0, the constant term is stored, at index 2 the coefficient
 * of the quadratic term is stored.
 */
void quadraticCoefs(FP_TYPE factor, int nReal, FP_TYPE a, FP_TYPE b, FP_TYPE *coefs) {
    switch (nReal) {
        case 0:
            // Two complex conjugate solutions a ± i*b.
            coefs[0] = factor * (a*a + b*b);
            coefs[1] = -2.0 * factor * a;
            coefs[2] = factor;
            break;
        case 2:
            // Two real solutions a and b.
            coefs[0] = factor * a * b;
            coefs[1] = -factor * (a + b);
            coefs[2] = factor;
            break;
        default:
            fprintf(stderr, "Real quadratic equation has either two real solutions or no real solution.\n");
            exit(1);
    }
}





// This computes (a*b - c*d) without intermediate rounding of the expressions a*b
// and c*d. Rounding only occurs at the end of the complete expression (a*b - c*d).
static inline FP_TYPE diff_of_products(FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE d) {
#if FP_PREC == 53
    FP_TYPE w = d * c;
    FP_TYPE e = fma(-d, c, w);
    FP_TYPE f = fma(a, b, -w);
    return f + e;
#endif
#if FP_PREC == 64
    return a*b - c*d;
#endif
}




// Like Java's Math.getExponent() method.
static inline int getExponent(FP_TYPE a) {
    int exp;
    FREXP(a, &exp);
    return exp - 1;
}

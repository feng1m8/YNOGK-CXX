#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "solve_utils.h"

static FP_TYPE diff_of_products(FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE d);
static FP_TYPE sum_of_products(FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE d);
static int getExponent(FP_TYPE a);


/**
 *
 * @author woelen
 */
#define PREC_THRESH  (-3*getExponent(EPS)/2)



/**
 * Solves a cubic equation with real coefficients of the form a*x^3 + b*x^2 + c*x + d == 0.
 * The complete solution is returned in the supplied arrays, in no particular order.
 * @param a The coefficient of the cubic term
 * @param b The coefficient of the quadratic term
 * @param c The coefficient of the linear term
 * @param d The constant term
 * @param zeror This is used for output. The first three elements of this array are
 * overwritten with the real part of the three zeros of the equation.
 * @param zeroi This is used for output. The first three elements of this array are
 * overwritten with the imaginary part of the three zeros of the equation.
 */
void cubicSolve(FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE d, FP_TYPE *zeror, FP_TYPE *zeroi) {

    if (d == 0.0) {
        quadraticSolve(a, b, c, zeror, zeroi);
        zeror[2] = zeroi[2] = 0.0;
        return;
    }

    FP_TYPE x = -b/(3*a);
    FP_TYPE q, dq;    // Polynomial value at x and its derivative at x
    FP_TYPE b1, c2;
    {
#if FP_PREC == 64
        FP_TYPE q0 = a*x;
        b1 = q0 + b;
        c2 = b1*x + c;
        dq = (q0 + b1)*x + c2;
        q = c2*x + d;
#endif
#if FP_PREC == 53
        b1 = fma(a, x, b);
        c2 = fma(b1, x, c);
        dq = fma(fma (a, x, b1), x, c2);
        q = fma(c2, x, d);
#endif
    }

    FP_TYPE r;
    int s;
    {
        FP_TYPE t = q/a;
        r = CBRT(abs(t));
        s = (t >= 0) ? 1 : -1;
        t = -dq/a;
        if (t > 0) {
            FP_TYPE L = 1.324717957244746;   // solution of L³ = L + 1
            r = L * max(r, SQRT(t));
        }
    }

    FP_TYPE x0 = x - s*r;   // The initial value for the Newton iteration.
    if (abs(x-x0) > EPS*abs(x)) {
        FP_TYPE SAFETY_FACTOR = 1.0 - 100*EPS;
        do {
            x = x0;
            {
#if FP_PREC == 64
                FP_TYPE q0 = a*x;
                b1 = q0 + b;
                c2 = b1*x + c;
                dq = (q0 + b1)*x + c2;
                q = c2*x + d;
#endif
#if FP_PREC == 53
                b1 = fma(a, x, b);
                c2 = fma(b1, x, c);
                dq = fma(fma (a, x, b1), x, c2);
                q = fma(c2, x, d);
#endif
            }
            x0 = (dq == 0) ? x : (x - ((q/dq)*SAFETY_FACTOR));
        } while ((s > 0) ? (x0 > x) : (x0 < x));

        FP_TYPE d_x = d/x;
        if (abs(a*x*x) > abs(d_x)) {
            c2 = -d_x;
            b1 = (c2 - c)/x;
        }
    }
    quadraticSolve(a, b1, c2, zeror, zeroi);
    zeror[2] = x;
    zeroi[2] = 0.0;
}





/**
 * Solves a cubic equation with real coefficients of the form x^3 + a*x^2 + b*x + c == 0.
 * The complete solution is returned in the supplied arrays, in no particular order.
 * @param a The coefficient of the quadratic term
 * @param b The coefficient of the linear term
 * @param c The constant term
 * @param zeror This is used for output. The first three elements of this array are
 * overwritten with the real part of the three zeros of the equation.
 * @param zeroi This is used for output. The first three elements of this array are
 * overwritten with the imaginary part of the three zeros of the equation.
 */
void cubicSolveMonic(FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE *zeror, FP_TYPE *zeroi) {

    if (c == 0.0) {
        quadraticSolveMonic(a, b, zeror, zeroi);
        zeror[2] = zeroi[2] = 0.0;
        return;
    }

    FP_TYPE x = -a/3;
    FP_TYPE q, dq;    // Polynomial value at x and its derivative at x
    FP_TYPE b1, c2;
    {         
        b1 = x + a;
#if FP_PREC == 64
        c2 = b1*x + b;
        dq = (x + b1)*x + c2;
        q = c2*x + c;
#endif
#if FP_PREC == 53
        c2 = fma(b1, x, b);
        dq = fma(x + b1, x, c2);
        q = fma(c2, x, c);
#endif
    }

    FP_TYPE r;
    int s;
    {
        r = CBRT(abs(q));
        s = (q >= 0) ? 1 : -1;
        if (dq < 0) {
            FP_TYPE L = 1.324717957244746;   // solution of L³ = L + 1
            r = L * max(r, SQRT(-dq));
        }
    }

    FP_TYPE x0 = x - s*r;   // The initial value for the Newton iteration.
    if (abs(x-x0) > EPS*abs(x)) {
        FP_TYPE SAFETY_FACTOR = 1.0 - 100*EPS;
        do {
            x = x0;
            {
                b1 = x + a;
#if FP_PREC == 64
                c2 = b1*x + b;
                dq = (x + b1)*x + c2;
                q = c2*x + c;
#endif
#if FP_PREC == 53
                c2 = fma(b1, x, b);
                dq = fma(x + b1, x, c2);
                q = fma(c2, x, c);
#endif
            }
            x0 = (dq == 0) ? x : (x - ((q/dq)*SAFETY_FACTOR));
        } while ((s > 0) ? (x0 > x) : (x0 < x));

        FP_TYPE d_x = c/x;
        if (x*x > abs(d_x)) {
            c2 = -d_x;
            b1 = (c2 - b)/x;
        }
    }
    quadraticSolveMonic(b1, c2, zeror, zeroi);
    zeror[2] = x;
    zeroi[2] = 0.0;
}




/**
 * Solves a cubic equation with real coefficients of the form x^3 + g*x + h == 0.
 * The complete solution is returned in the supplied arrays, in no particular order.
 * @param g The coefficient of the linear term
 * @param h The constant term
 * @param zeror This is used for output. The first three elements of this array are
 * overwritten with the real part of the three zeros of the equation.
 * @param zeroi This is used for output. The first three elements of this array are
 * overwritten with the imaginary part of the three zeros of the equation.
 */
void cubicSolveDepressed(FP_TYPE g, FP_TYPE h, FP_TYPE *zeror, FP_TYPE *zeroi) {
    FP_TYPE x = cubicSolveDepressedDominant(g, h);
    FP_TYPE xx = x*x;
    FP_TYPE b = (g >= 0 || g < -2*xx) ? g + xx : -h/x;
    quadraticSolveMonic(x, b, zeror, zeroi);
    zeror[2] = x;
    zeroi[2] = 0.0;
}




//*******************************************************************
//  solveDepressedDominant: Finds the dominant root of a depressed cubic
//               q(x) = x^3 + g*x + h with real coefficients g and h,
//               using the parabola/reciprocal intersection method,
//               introduced by Peter Strobach.
//               The algorithm is enhanced with scaling to make it
//               robust against underflow and overflow.
//====================================================================


/**
 * Solves an equation of the form x^3 + g*x + h == 0. The so-called
 * dominant solution is returned. This is the real solution if there
 * only is one real solution, and it is the best conditioned real
 * solution if there are three solutions. This method is intended to
 * be used for special purposes, and it also is used for solving
 * quartic equations. This method only returns one solution of the
 * depressed equation. If you want the full solution, then use the
 * other solveDepressed(...) method, which returns the full set of
 * all three zeros.
 * @param g The coefficient of the linear term
 * @param h The constant term
 * @return The (real) dominant solution of the depressed cubic equation.
 */
FP_TYPE cubicSolveDepressedDominant(FP_TYPE g, FP_TYPE h) {
    // First check some extreme situations, which can be handled
    // easily, and also can be handled in such a way that no
    // overflow or underflow occurs.
    int expH = getExponent(h)/3;
    int dExpG = getExponent(g) - 2*expH;
    if (g == 0.0 || dExpG < -PREC_THRESH) {
        return CBRT(-h);
    }
    if (h == 0.0 || dExpG > PREC_THRESH) {
        return (g > 0) ? (-h/g) : SQRT(-g);
    }

    // Now we know that the coefficients can safely be scaled
    // in order to avoid underflow or overflow. Scale them, then
    // solve the equation, and at the end correct the solution.
    g = SCALB(g, -2*expH);
    h = SCALB(h, -3*expH);

    // Here, the original algorithm from Peter Strobach is started.

    FP_TYPE x0, xr;
    if (h < 0.0) {
        xr = SQRT(-h);               // equation (3.20) of Strobachs paper
        if (g > xr) {
            x0 = -h / g;             // equation (3.26) of Strobachs paper
        } 
        else if (g < h) {            // g < -xr * xr
            x0 = SQRT(-g);           // equation (3.22) of Strobachs paper
        } 
        else {
            x0 = xr;                 // equation (3.24) of Strobachs paper
        }
    }
    else {    // (h > 0.0)
        xr = SQRT(h);                // equation (3.21) of Strobachs paper
        if (g > xr) {
            x0 = -h / g;             // equation (3.26) of Strobachs paper
        } 
        else if (g < -h) {           // g < -xr * xr
            x0 = -SQRT(-g);          // equation (3.23) of Strobachs paper
        } 
        else {
            x0 = -xr;                // equation (3.25) of Strobachs paper
        }
    }

    // Perform a set of so-called "bi-iterations", as described
    // in Strobachs paper, just after equation (3.11). This consists
    // of two steps, both are programmed below.
    FP_TYPE xDiff = MAX_VALUE;
    for (int iter = 0; x0==x0 && iter < 8; iter++) {
        FP_TYPE a, b, c, xmin, xmax;

        // Compute coefficients of equation (3.10) of Strobachs paper
        a = x0 * x0;
        b = -h;
        c = sum_of_products(g, a, 2*h, x0);

        // Solve Strobachs equation (3.10), determine the smallest
        // solution and put that in the next iteration of x0.
        if (a == 0.0) {
            xmin = 0.0;
        }
        else {
            int exp = -max(getExponent(b), (2 + getExponent(a) + getExponent(c))/2);
            FP_TYPE D = diff_of_products(SCALB(b, exp), SCALB(b, exp), 4*SCALB(a, exp), SCALB(c, exp));
            if (D <= 0.0) {
                xmin = -b/(2*a);  // Only the real part is used.
            }
            else {
                FP_TYPE num = (b > 0.0) ? (-b - SCALB(SQRT(D), -exp)) : (-b + SCALB(SQRT(D), -exp));
                xmin = (2*c) / num;  // smallest magnitude also is lowest value.
            }    
        }
        x0 = xmin;  // equation (3.11) of Strobachs paper

        // Here we are at the end of step 1 of the "bi-iteration"
        // and we perform the second step of the " bi-iteration".


        // Compute coefficients of equation (3.15) of Strobachs paper
        a = 2 * x0;
#if FP_PREC == 64
        b = -x0*x0 + g;
#endif
#if FP_PREC == 53
        b = fma(x0, -x0, g);
#endif
        c = h;

        if (a == 0.0) {
            xmin = xmax = 0.0;
        }
        else {
            int exp = -max(getExponent(b), (2 + getExponent(a) + getExponent(c))/2);
            FP_TYPE D = diff_of_products(SCALB(b, exp), SCALB(b, exp), 4*SCALB(a, exp), SCALB(c, exp));
            if (D <= 0.0) {
                xmin = xmax = -b/(2*a);
            }
            else {
                FP_TYPE num = (b > 0.0) ? (-b - SCALB(SQRT(D), -exp)) : (-b + SCALB(SQRT(D), -exp));
                xmax = num / (2*a);   // highest magnitude, not necessarily highest value
                xmin = (2*c) / num;   // lowest magnitude, not necessarily lowest value
            }    
        }

        FP_TYPE x0Prev = x0;

        // At this point we have two real solutions. One of
        // these is positive and one of these is negative.
        // Use equations (3.16) and (3.17) of Strobachs paper
        // to determine the new value of x0.
        x0 = (h < 0)
                  ? max(xmin, xmax)     // take the positive solution: eq. (3.16) of Strobachs paper
                  : min(xmin, xmax);    // take the negative solution: eq. (3.17) of Strobachs paper

        FP_TYPE xDiffPrev = xDiff;
        xDiff = abs(x0Prev - x0);
        if (xDiff == 0.0 || xDiff == xDiffPrev) {
            // terminate
            break;
        }
    }

    // Return the solution x0, but first undo the effect of scaling.
    return SCALB(x0, expH);
}


    
    
    
    
    
/**
 * Given a set of roots and the coefficient of the leading term, compute
 * all coefficients of the cubic equation.
 * @param factor The coefficient of the leading term (the cubic term).
 * @param nReal Tells whether the equation has three real roots, or one real
 * root. The next three parameters have meaning, depending on this parameter.
 * @param a This is interpreted as a real root, regardless of the value of
 * nReal.
 * @param b If nReal equals 3, then this is interpreted as one of the real
 * roots of the equation. If nReal equals 1, then this is interpreted as the
 * real part of both complex conjugate roots.
 * @param c If nReal equals 3, then this is interpreted as one of the real
 * roots of the equation. If nReal equals 1, then this is interpreted as the
 * imaginary part of one of the complex conjugate roots, the other complex
 * root has imaginary part equal to -c.
 * @param coefs An output array, which is filled with the coefficients
 * of the cubic. This array must have a length at least equal to 4. 
 * At index 0, the constant term is stored, at index 3 the coefficient
 * of the cubic term is stored.
 */
void cubicCoefs(FP_TYPE factor, int nReal, FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE *coefs) {
    // First compute the quadratic factor x^2 + p*x + q, based on b and c.
    FP_TYPE p, q;
    switch (nReal) {
        case 1:
            // One real solution a and two complex conjugate solutions b ± i*c.
            p = -2.0 * b;
            q = b*b + c*c;
            break;
        case 3:
            // Three real solutions a, b, and c. The arguments b and c are
            // regarded as real roots of the equation.
            p = -(b + c);
            q = b * c;
            break;
        default:
            fprintf(stderr, "Real cubic equation has either three real solutions or one real solution.\n");
            exit(1);
    }

    // Compute factor * (x - a) * (x^2 + px + q).
    coefs[0] = -factor * a * q;
    coefs[1] = factor * (q - a*p);
    coefs[2] = factor * (p - a);
    coefs[3] = factor;
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







// This computes (a*b + c*d) without intermediate rounding of the expressions a*b
// and c*d. Rounding only occurs at the end of the complete expression (a*b + c*d).
static inline FP_TYPE sum_of_products(FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE d) {
#if FP_PREC == 53
    FP_TYPE w = d * c;
    FP_TYPE e = fma(d, c, -w);
    FP_TYPE f = fma(a, b, w);
    return f + e;
#endif
#if FP_PREC == 64
    return a*b + c*d;
#endif
}




// Like Java's Math.getExponent() method.
static inline int getExponent(FP_TYPE a) {
    int exp;
    FREXP(a, &exp);
    return exp - 1;
}





/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

// Code below is left in for documentation purposes. The code above is
// optimized and made more accurate, using FMA and scaling. The code
// below is the basic algorithm, presented by Peter Strobach.

// Find the dominant real root of the cubic x^3 + g*x + h,
// using the parabola/reciprocal intersection method, introduced
// by Peter Strobach.
/*
static FP_TYPE depressed_cubic_root(FP_TYPE g, FP_TYPE h) {
    if (h == 0.0) {
        return (g < 0) ? SQRT(-g) : 0.0;
    } 

    FP_TYPE x0, xr;
    if (h < 0.0) {
        xr = SQRT(-h);   // equation (3.20) of Strobachs paper
        if (g > xr) {
            x0 = -h / g;             // equation (3.26) of Strobachs paper
        } 
        else if (g < h) {            // g < -xr * xr
            x0 = SQRT(-g);      // equation (3.22) of Strobachs paper
        } 
        else {
            x0 = xr;                 // equation (3.24) of Strobachs paper
        }
    }
    else {  
        // (h > 0.0)
        xr = SQRT(h);    // equation (3.21) of Strobachs paper
        if (g > xr) {
            x0 = -h / g;             // equation (3.26) of Strobachs paper
        } 
        else if (g < -h) {           // g < -xr * xr
            x0 = -SQRT(-g);     // equation (3.23) of Strobachs paper
        } 
        else {
            x0 = -xr;                // equation (3.25) of Strobachs paper
        }
    }

    // Perform a set of so-called "bi-iterations", as described
    // in Strobachs paper, just after equation (3.11). This consists
    // of two steps, both are programmed below.
    FP_TYPE xDiff = MAX_VALUE;
    for (int iter = 0; x0==x0 && iter < 8; iter++) {
        FP_TYPE a, b, c, xmin, xmax;

        // Compute coefficients of equation (3.10) of Strobachs paper
        a = x0 * x0;
        b = -h;
        c = g * a + 2 * h * x0;

        // Solve Strobachs equation (3.10), determine the smallest
        // solution and put that in the next iteration of x0.
        if (a == 0.0) {
            xmin = 0.0;
        }
        else {
            FP_TYPE D = b*b - 4*a*c;
            if (D < 0.0) {
                xmin = -b/(2*a);  // Only the real part is used.
            }
            else {
                FP_TYPE num = (b > 0.0) ? (-b - SQRT(D)) : (-b + SQRT(D));
                xmin = (2*c) / num;  // smallest magnitude also is lowest value.
            }    
        }
        x0 = xmin;  // equation (3.11) of Strobachs paper

        // Here we are at the end of step 1 of the "bi-iteration"
        // and we perform the second step of the " bi-iteration".


        // Compute coefficients of equation (3.15) of Strobachs paper
        a = 2 * x0;
        b = g - x0 * x0;
        c = h;

        if (a == 0.0) {
            xmin = xmax = 0.0;
        }
        else {
            FP_TYPE D = b*b - 4*a*c;
            if (D < 0.0) {
                xmin = xmax = -b/(2*a);
            }
            else {
                FP_TYPE num = (b > 0.0) ? (-b - SQRT(D)) : (-b + SQRT(D));
                xmax = num / (2*a);   // highest magnitude, not necessarily highest value
                xmin = (2*c) / num;   // lowest magnitude, not necessarily lowest value
            }    
        }

        FP_TYPE x0Prev = x0;

        // At this point we have two real solutions. One of
        // these is positive and one of these is negative.
        // Use equations (3.16) and (3.17) of Strobachs paper
        // to determine the new value of x0.
        x0 = (h < 0)
                  ? max(xmin, xmax)     // take the positive solution: eq. (3.16) of Strobachs paper
                  : min(xmin, xmax);    // take the negative solution: eq. (3.17) of Strobachs paper

        FP_TYPE xDiffPrev = xDiff;
        xDiff = abs(x0Prev - x0);
        if (xDiff == 0.0 || xDiff == xDiffPrev) {
            // terminate
            break;
        }
    }

    return x0;
}
*/
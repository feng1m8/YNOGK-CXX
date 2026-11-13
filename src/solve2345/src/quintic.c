#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "solve_utils.h"

#define PREC_THRESH  (-5*getExponent(EPS)/4)


static void roots(FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE d, FP_TYPE e, FP_TYPE f,
                  FP_TYPE *rootsRe, FP_TYPE *rootsIm);
static FP_TYPE bisect(FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE d, FP_TYPE e, FP_TYPE f);

static void cdiv(FP_TYPE numRe, FP_TYPE numIm, FP_TYPE denRe, FP_TYPE denIm, FP_TYPE *resRe, FP_TYPE *resIm);
static FP_TYPE dhornerReal(FP_TYPE da, FP_TYPE db, FP_TYPE dc, FP_TYPE dd, FP_TYPE de, FP_TYPE z);
static void dhorner(FP_TYPE da, FP_TYPE db, FP_TYPE dc, FP_TYPE dd, FP_TYPE de, FP_TYPE zRe, FP_TYPE zIm, FP_TYPE *resRe, FP_TYPE *resIm);
static FP_TYPE hornerReal(FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE d, FP_TYPE e, FP_TYPE f, FP_TYPE z);
static void horner(FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE d, FP_TYPE e, FP_TYPE f, FP_TYPE zRe, FP_TYPE zIm, FP_TYPE *resRe, FP_TYPE *resIm);
static FP_TYPE closestRootDistance(int idx, FP_TYPE *rootsRe, FP_TYPE *rootsIm);
static int getExponent(FP_TYPE a);

#if FP_PREC == 64
static inline FP_TYPE hypotl_impl(FP_TYPE a, FP_TYPE b) {return SQRT(a*a+b*b);} 
#endif

/**
 *
 * @author woelen
 */


/**
 * Solves a quintic equation with real coefficients of the form a*x^5 + b*x^4 + c*x^3 + d*x^2 + e*x + f == 0.
 * The complete solution is returned in the supplied arrays, in no particular order.
 * @param a The coefficient of the quintic term
 * @param b The coefficient of the quartic term
 * @param c The coefficient of the cubic term
 * @param d The coefficient of the quadratic term
 * @param e The coefficient of the linear term
 * @param f The constant term
 * @param zeror This is used for output. The first five elements of this array are
 * overwritten with the real part of the five zeros of the equation.
 * @param zeroi This is used for output. The first five elements of this array are
 * overwritten with the imaginary part of the five zeros of the equation.
 */    
void quinticSolve(FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE d, FP_TYPE e, FP_TYPE f,
                  FP_TYPE *zeror, FP_TYPE *zeroi) {
    if (f == 0) {
        zeror[4] = zeroi[4] = 0.0;
        quarticSolve(a, b, c, d, e, zeror, zeroi);
        return;
    }

    roots(a, b, c, d, e, f, zeror, zeroi);
}






/**
 * Solves a quintic equation with real coefficients of the form 
 * x^5 + a*x^4 + b*x^3 + c*x^2 + d*x + e == 0.
 * The complete solution is returned in the supplied arrays, in no particular order.
 * @param a The coefficient of the quartic term
 * @param b The coefficient of the cubic term
 * @param c The coefficient of the quadratic term
 * @param d The coefficient of the linear term
 * @param e The constant term
 * @param zeror This is used for output. The first five elements of this array are
 * overwritten with the real part of the five zeros of the equation.
 * @param zeroi This is used for output. The first five elements of this array are
 * overwritten with the imaginary part of the five zeros of the equation.
 */  
void quinticSolveMonic(FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE d, FP_TYPE e,
                       FP_TYPE *zeror, FP_TYPE *zeroi) {
    if (e == 0) {
        zeror[4] = zeroi[4] = 0.0;
        quarticSolveMonic(a, b, c, d, zeror, zeroi);
        return;
    }

    roots(1.0, a, b, c, d, e, zeror, zeroi);
}






/**
 * Solves an equation of the form x^5 + g*x + h == 0. The so-called
 * dominant solution is returned. This is the real solution if there
 * only is one real solution, and it is the best conditioned real
 * solution if there are three real solutions. This method is intended
 * to be used for special purposes. This method only returns one solution
 * of the depressed equation. If you want the full solution, then use
 * the other solveBringJerrard(...) method, which returns the full
 * set of all five zeros.
 * @param g The coefficient of the linear term
 * @param h The constant term
 * @return The (real) dominant solution of the fully depressed quintic
 * equation.
 */
FP_TYPE quinticSolveBringJerrardDominant(FP_TYPE g, FP_TYPE h) {
    // First check some extreme situations, which can be handled
    // easily, and also can be handled in such a way that no
    // overflow or underflow occurs.
    int expH = getExponent(h)/5;
    int dExpG = getExponent(g) - 4*expH;
    if (g == 0.0 || dExpG < -PREC_THRESH) {
        FP_TYPE x = (h > 0) ? -POW(h, 0.2) : POW(-h, 0.2);
        // We perform one Newton step (with g left out), in order
        // to improve the quality of the root. The computation of
        // pow() is somewhat inaccurate (2 ulps or so).
        FP_TYPE xx = x*x;
        FP_TYPE xxxx = xx*xx;
        FP_TYPE corr = (x*xxxx+h)/(5*xxxx);
        return x - corr;

    }
    if (h == 0.0 || dExpG > PREC_THRESH) {
        return (g > 0) ? -h/g : SQRT(SQRT(-g));
    }

    // Now we know that the coefficients can safely be scaled
    // in order to avoid underflow or overflow. Scale them, then
    // solve the equation, and at the end correct the solution.
    g = SCALB(g, -4*expH);
    h = SCALB(h, -5*expH);

    // Determine an initial estimate of x, used for the Newton iteration.
    // Use an improved Lagrange upper bound for the magnitude of a root,
    // and for g > 0, also take into account |h/g| as a bound.
    FP_TYPE x = POW(abs(h), 0.2) + SQRT(SQRT(abs(g)));
    if (g > 0) {
        x = min(x, abs(h)/g);
    }
    if (h > 0) {
        x = -x;
    }

    // Perform the Newton iteration, until the correction
    // does not change x anymore, or until the sign of the
    // correction changes. The latter is a guard against
    // x flipping around the zero with very small steps.
    boolean corrSign = (h > 0);
    while (true) {
        FP_TYPE xx = x*x;
        FP_TYPE xxxx = xx*xx;
        FP_TYPE corr = (x*(xxxx + g) + h) / (5*xxxx + g);
        FP_TYPE xcorr = x - corr;
        if (xcorr == x || ((corr < 0) != corrSign)) {
            break;
        }
        x = xcorr;
    }

    // Return the solution x, but first undo the effect of scaling.
    return SCALB(x, expH);
}





/**
 * Solves a quintic equation with real coefficients of the form x^5 + g*x + h == 0.
 * The complete solution is returned in the supplied arrays, in no particular order.
 * @param g The coefficient of the linear term
 * @param h The constant term
 * @param zeror This is used for output. The first five elements of this array are
 * overwritten with the real part of the five zeros of the quation.
 * @param zeroi This is used for output. The first five elements of this array are
 * overwritten with the imaginary part of the five zeros of the quation.
 */
void quinticSolveBringJerrard(FP_TYPE g, FP_TYPE h, FP_TYPE *zeror, FP_TYPE *zeroi) {
    FP_TYPE x = quinticSolveBringJerrardDominant(g, h);
    FP_TYPE xx = x*x;
    FP_TYPE xxx = xx*x;
    FP_TYPE xxxx = xx*xx;
    FP_TYPE d = (g >= 0 || g < -2*xxxx) ? g + xxxx : -h/x;
    quarticSolveMonic(x, xx, xxx, d, zeror, zeroi);
    zeror[4] = x;
    zeroi[4] = 0.0;
}





    
    
    
    
    
    
/**
 * Given a set of roots and the coefficient of the leading term, compute
 * all coefficients of the quintic equation.
 * @param factor The coefficient of the leading term (the quintic term).
 * @param nReal Tells whether the equation has 1, 3, or 5 real roots. The 
 * next five parameters have meaning, depending on this parameter. The
 * first nReal parameters are regarded as real roots, the rest of these
 * five parameters are regarded as the real and imaginary part of complex
 * conjugate roots.
 * @param a This is interpreted as a real root of the quintic equation,
 * regardless of the value of nReal.
 * @param b If nReal equals 3 or 5, then this is interpreted as one of the
 * real roots of the equation. If nReal equals 1, then this is interpreted
 * as the real part of a pair of complex conjugate roots.
 * @param c If nReal equals 3 or 5, then this is interpreted as one of the
 * real roots of the equation. If nReal equals 1, then this is interpreted
 * as the imaginary part of a pair of complex conjugate roots, the other
 * conjugate complex root has imaginary part equal to -c.
 * @param d If nReal equals 5, then this is interpreted as one of the real
 * roots of the equation. If nReal equals 1 or 3, then this is interpreted
 * as the real part of a pair of complex conjugate roots.
 * @param e If nReal equals 5, then this is interpreted as one of the real
 * roots of the equation. If nReal equals 1 or 3, then this is interpreted
 * as the imaginary part of a pair of complex conjugate roots, the other
 * conjugate complex root has imaginary part equal to -e.
 * @param coefs An output array, which is filled with the coefficients
 * of the quintic. This array must have a length at least equal to 6. 
 * At index 0, the constant term is stored, at index 5 the coefficient
 * of the quintic term is stored.
 */
void quinticCoefs(FP_TYPE factor, int nReal, FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE d, FP_TYPE e, FP_TYPE *coefs) {
    // First, compute the quadratic factor x^2 + p*x + q, based on b and c.
    FP_TYPE p, q;
    switch (nReal) {
        case 1:
            // There is one real solution, the arguments b and c are regarded
            // as real and imaginary parts of the complex roots b ± i*c.
            p = -2.0 * b;
            q = b*b + c*c;
            break;
        case 3:
        case 5:
            // Either 3, or 5 real roots. The parameters b and c are
            // regarded as real roots.
            p = -(b + c);
            q = b * c;
            break;
        default:
            fprintf(stderr, "Real quintic equation has either five real solutions, three real solutions or one real solution.\n");
            exit(1);
    }

    // Next, compute the quadratic factor, based on d and e.
    FP_TYPE r, s;
    switch (nReal) {
        case 1:
        case 3:
            // There is only one real solution or three real solutions.
            // The arguments d and e are regarded as real and imaginary
            // parts of the complex roots d ± i*e.
            r = -2.0 * d;
            s = d*d + e*e;
            break;
        case 5:
            // There are 5 real roots. The parameters d and e are
            // regarded as real roots.
            r = -(d + e);
            s = d * e;
            break;
        default:
            fprintf(stderr, "Real quintic equation has either five real solutions, three real solutions or one real solution.\n");
            exit(1);
    }

    // Compute x^4 + k*x^3 + l*x^2 + m*x + n = (x^2 + px + q) * (x^2 + rx + s).
    FP_TYPE n = q * s;
    FP_TYPE m = p*s + q*r;
    FP_TYPE l = p*r + q + s;
    FP_TYPE k = p + r;

    // Compute factor * (x - a) * (x^4 + k*x^3 + l*x^2 + m*x + n).
    coefs[0] = -factor * a * n;
    coefs[1] = factor * (n - a*m);
    coefs[2] = factor * (m - a*l);
    coefs[3] = factor * (l - a*k);
    coefs[4] = factor * (k - a);
    coefs[5] = factor;
}





// This is the code, which does the actual root finding. It is called
// from within the wrapper method solve() or solveMonic().
//
// rootsRe, rootsIm: A pair of arrays, which must have a length
//    of at least 5 elements. These arrays will be filled by
//    the code. Any previous contents is overwritten.
//
static void roots(FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE d, FP_TYPE e, FP_TYPE f,
                  FP_TYPE *rootsRe, FP_TYPE *rootsIm) {

    int MAX_ITER = 100;

    FP_TYPE absf = abs(f);
    FP_TYPE abse = abs(e);
    FP_TYPE absd = abs(d);
    FP_TYPE absc = abs(c);
    FP_TYPE absb = abs(b);
    FP_TYPE absa = abs(a);

    // Scale coefficients
    {
        FP_TYPE u1 = max(absa, max(absb, max(absc, max(absd, max(abse, absf)))));
        FP_TYPE u2 = min(absa, absf);
        int exp1 = getExponent(u1);
        int exp2 = getExponent(u2);

        int exp = -(exp1 + exp2)/2;

        FP_TYPE u = SCALB(1.0, exp);
        a *= u;
        b *= u;
        c *= u;
        d *= u;
        e *= u;
        f *= u;
        absa *= u;
        absb *= u;
        absc *= u;
        absd *= u;
        abse *= u;
        absf *= u;
    }

    // Compute coefficients of derivative of polynomial
    FP_TYPE da = 5*a;
    FP_TYPE db = 4*b;
    FP_TYPE dc = 3*c;
    FP_TYPE dd = 2*d;
    FP_TYPE de = e;

    // The variables below are used for the following purposes:
    // zRe, zIm: Current iterate of the zero to be found.
    // z0Re, z0Im: Previous iterate of the zero to be found.
    // ztRe, ztIm: Tentative iterate of the zero to be found.
    // wRe, wIm: Current search point along a line for finding
    //    the best next iterate (with smallest value of f).
    // afz: |f(Z)|, where Z=zRe+i*zIm is the current iterate of the
    //    zero, and f stands for the polynomial.
    // afzmin: Estimate of the lowest attainable value of afz.
    // afz0: |f(Z0)|, where Z0=z0Re+i*z0Im is the previous iterate
    //    of the zero.
    // dzRe, dzIm: Step to be taken in zRe and zIm. This is a tentative
    //    step, which can be halved, FP_TYPEd, or incremented in equal
    //    steps along a line (see paper PDF file).
    // fzRe, fzIm: f(Z), real part and imaginary part.
    // f1zRe, f1zIm: f'(Z), real part and imaginary part.
    // f1z0Re, f1z0Im: f'(Z0), real part and imaginary part.
    // nn: Order of the deflated polynomial.
    // stage1: Tells whether for the current deflated polynomial the
    //    algorithm is in the searching stage (along a line wRe, wIm)
    //    or in a mode, where simple Newton iteration gives a good
    //    and certain next iterate. When in searching mode, this is
    //    equal to true.
    FP_TYPE wRe, wIm, zRe, zIm, z0Re, z0Im, ztRe, ztIm;
    FP_TYPE dzRe, dzIm, f1zRe, f1zIm, f1z0Re, f1z0Im, fzRe, fzIm;
    FP_TYPE afz, afz0, afzmin;
    FP_TYPE r0;
    int iter;
    boolean stage1 = false;

    // Set initial iterates and bound used in the convergence test.
    z0Re = z0Im = 0.0;
    afz0 = absf;
    afzmin = afz0 * 5 * 16.0 * EPS;

    ztRe = z0Re;
    ztIm = z0Im;
    r0 = POW(absf/absa, 0.2);
    if (b != 0.0) r0 = min(r0, SQRT(SQRT(absf/absb)));
    if (c != 0.0) r0 = min(r0, CBRT(absf/absc));
    if (d != 0.0) r0 = min(r0, SQRT(absf/absd));
    if (e != 0.0) r0 = min(r0, absf/abse);
    r0 = r0/2;

    f1z0Re = e;
    f1z0Im = 0.0;
    zRe = ((f > 0.0) == (e > 0.0)) ? -r0 : r0; 
    zIm = 0.0;

    dzRe = zRe;
    dzIm = zIm;
    fzRe = hornerReal(a, b, c, d, e, f, zRe);
    fzIm = 0.0;
    afz = abs(fzRe);

    for (iter = 0; iter < MAX_ITER; iter++) {
        if (iter == 0 || afz < afz0) {
            // First iteration or iteration following a successful step.
            // Calculate the tentative step DZ and whether in stage 1.
            dhorner(da, db, dc, dd, de, zRe, zIm, &f1zRe, &f1zIm);
            FP_TYPE af1z = HYPOT(f1zRe, f1zIm);
            if (af1z != 0.0) {
                cdiv(fzRe, fzIm, -f1zRe, -f1zIm, &dzRe, &dzIm);
                FP_TYPE af2z = HYPOT(f1z0Re-f1zRe, f1z0Im-f1zIm)/HYPOT(z0Re-zRe, z0Im-zIm);  // Appr. abs(f''(z))
                stage1 = afz*af2z/af1z > af1z*0.5 || zRe != ztRe || zIm != ztIm;
                FP_TYPE r = HYPOT(dzRe, dzIm);
                if (r > r0 * 3.0) {
                    FP_TYPE factor = r0/r;
                    FP_TYPE re = dzRe*1.8 - dzIm*2.4;
                    FP_TYPE im = dzRe*2.4 + dzIm*1.8;
                    dzRe = re*factor;
                    dzIm = im*factor;
                }
            } 
            else {
                FP_TYPE re = dzRe*1.8 - dzIm*2.4;
                FP_TYPE im = dzRe*2.4 + dzIm*1.8;
                dzRe = re;
                dzIm = im;
                stage1 = true;
            }

            f1z0Re = f1zRe;
            f1z0Im = f1zIm;
        }

        // Find the next point in the iteration.
        z0Re = zRe;
        z0Im = zIm;
        afz0 = afz;
        zRe = z0Re + dzRe;
        zIm = z0Im + dzIm;
        wRe = zRe;
        wIm = zIm;
        horner(a, b, c, d, e, f, zRe, zIm, &fzRe, &fzIm);
        afz = HYPOT(fzRe, fzIm);
        ztRe = zRe;
        ztIm = zIm;
        if (stage1) {
            // Beginning of stage 1 search.
            // If div2 is true, then we search with steps w = dz, dz/2, dz/4, ...
            // If div2 is false, then we search with steps w = dz, 2*dz, 4*dz, ...
            boolean div2 = (afz >= afz0);
            for (int k = 1; k <= 5; k++) {
                if (div2) {
                    if (k <= 2) {
                        wRe = (wRe + z0Re)*0.5;
                        wIm = (wIm + z0Im)*0.5;
                    } 
                    else if (k == 3) {
                        FP_TYPE re = dzRe*0.15 - dzIm*0.20;
                        FP_TYPE im = dzRe*0.20 + dzIm*0.15;
                        wRe = z0Re + re;
                        wIm = z0Im + im;
                    } 
                    else {
                        break;
                    }
                } 
                else {
                    wRe += dzRe;
                    wIm += dzIm;
                }

                FP_TYPE fwRe, fwIm;
                horner(a, b, c, d, e, f, wRe, wIm, &fwRe, &fwIm);
                FP_TYPE afw = HYPOT(fwRe, fwIm);
                if (afw >= afz) {
                    break;
                }
                afz = afw;
                fzRe = fwRe;
                fzIm = fwIm;
                zRe = wRe;
                zIm = wIm;
            } // end of k-loop
        }
        r0 = HYPOT(z0Re-zRe, z0Im-zIm);

        // Test whether we converged to a zero or not.
        if (afz >= afz0) {
            zRe = z0Re;
            zIm = z0Im;
        }
        if (r0 < EPS * HYPOT(zRe, zIm)) {
            FP_TYPE afaz = hornerReal(absa, absb, absc, absd, abse, absf, HYPOT(zRe, zIm));
            if (afz < 16*5*EPS*afaz) {
		// We converged, we found a root.
                break;
            }

            // We're in a very flat area around the origin,
            // escape from that area carefully.
            do {
                zRe *= 1.1;
                zIm *= 1.1;
                horner(a, b, c, d, e, f, zRe, zIm, &fzRe, &fzIm);
                afz0 = afz;
                afz = HYPOT(fzRe, fzIm);
            }
            while (abs(afz-afz0) < 0.1*(afz + afz0));

            // We're out, and change the search direction.
            dzRe = -0.1*zIm;
            dzIm = 0.1*zRe;
            z0Re = zRe - dzRe;
            z0Im = zIm - dzIm;
            continue;
        }
        if (afz*1.0001 >= afz0) {
            afz = afz0;
            if (afz <= afzmin) {
                break;
            }
            // We did not converge successfully. Make the tentative step
            // half as long and change its direction in the complex plane.
            FP_TYPE re = dzRe*0.3 - dzIm*0.4;
            FP_TYPE im = dzRe*0.4 + dzIm*0.3;
            dzRe = -re;
            dzIm = -im;
        }
    }  // end of iter loop

    if (iter == MAX_ITER || zRe != zRe) {
        // No convergence. We have a quintic equation with real coefficients,
        // so we can always use interval bisection to find a real root.
        // It is relatively inefficient, but it is the best we can do.
        // At least we do find roots  and do not have to give up.
        zRe = bisect(a, b, c, d, e, f);
        zIm = 0;
    }

    // Deflate polynomial coefficients, and store the root.
    // If the root is a conjugate complex root, then store
    // two roots and deflate with a quadratic factor.
    if (abs(zIm) < 5*EPS*abs(zRe)) {
        // Deflation with a single real root, 100% safe to do so.
        rootsRe[4] = zRe;
        rootsIm[4] = 0.0;

        // Always perform forward deflation. This works fine, because
        // roots are roughly found in order of increasing magnitude.
#if FP_PREC == 53
        FP_TYPE a4 = a;
        FP_TYPE b4 = fma(a4, zRe, b);
        FP_TYPE c4 = fma(b4, zRe, c);
        FP_TYPE d4 = fma(c4, zRe, d);
        FP_TYPE e4 = fma(d4, zRe, e);
#endif
#if FP_PREC == 64
        FP_TYPE a4 = a;
        FP_TYPE b4 = a4*zRe + b;
        FP_TYPE c4 = b4*zRe + c;
        FP_TYPE d4 = c4*zRe + d;
        FP_TYPE e4 = d4*zRe + e;
#endif
        quarticSolve(a4, b4, c4, d4, e4, rootsRe, rootsIm);
    } 
    else {
        // We found a root, which may be one of a conjugate pair,
        // but it also may be a real root, into which a small
        // imaginary part has crept during the iteration. Only
        // if we have great confidence that this root is one of
        // a conjugate pair, we deflate with a quadratic factor,
        // otherwise we deflate with a real root, accepting a
        // possible small error, which may affect another root.
        FP_TYPE abszRe = abs(zRe);
        FP_TYPE g = 0.0;
        FP_TYPE fc = a;
#if FP_PREC == 53
        {
            FP_TYPE fd = fma(fc, zRe, b);
            g = abszRe * (g + abs(fc)) + abs(fd);
            fc = fd;
        }
        {
            FP_TYPE fd = fma(fc, zRe, c);
            g = abszRe * (g + abs(fc)) + abs(fd);
            fc = fd;
        }
        {
            FP_TYPE fd = fma(fc, zRe, d);
            g = abszRe * (g + abs(fc)) + abs(fd);
            fc = fd;
        }
        {
            FP_TYPE fd = fma(fc, zRe, e);
            g = abszRe * (g + abs(fc)) + abs(fd);
            fc = fd;
        }
        {
            FP_TYPE fd = fma(fc, zRe, f);
            g = abszRe * (g + abs(fc)) + abs(fd);
            fc = fd;
        }
#endif
#if FP_PREC == 64
        {
            FP_TYPE fd = fc*zRe + b;
            g = abszRe * (g + abs(fc)) + abs(fd);
            fc = fd;
        }
        {
            FP_TYPE fd = fc*zRe + c;
            g = abszRe * (g + abs(fc)) + abs(fd);
            fc = fd;
        }
        {
            FP_TYPE fd = fc*zRe + d;
            g = abszRe * (g + abs(fc)) + abs(fd);
            fc = fd;
        }
        {
            FP_TYPE fd = fc*zRe + e;
            g = abszRe * (g + abs(fc)) + abs(fd);
            fc = fd;
        }
        {
            FP_TYPE fd = fc*zRe + f;
            g = abszRe * (g + abs(fc)) + abs(fd);
            fc = fd;
        }
#endif
        if (abs(fc) > 2.0 * EPS * g + min(afz, afz0)) {
            // We have great confidence that the root is one of
            // a pair of conjugate roots. We deflate with a 
            // quadratic factor and add two roots to the output.
            rootsRe[3] = rootsRe[4] = zRe;
            rootsIm[3] = -(rootsIm[4] = zIm);

            FP_TYPE pp = -zRe - zRe;
            FP_TYPE qq = zRe * zRe + zIm * zIm;
            FP_TYPE a3 = a;
            FP_TYPE b3 = b - pp * a3;
            FP_TYPE c3 = c - pp * b3 - qq * a3;
            FP_TYPE d3 = d - pp * c3 - qq * b3;
            cubicSolve(a3, b3, c3, d3, rootsRe, rootsIm);
        } 
        else {
            // We deflate with a single root and only
            // use the real part of the root we found.
            rootsRe[4] = zRe;
            rootsIm[4] = 0.0;
#if FP_PREC == 53
            FP_TYPE a4 = a;
            FP_TYPE b4 = fma(a4, zRe, b);
            FP_TYPE c4 = fma(b4, zRe, c);
            FP_TYPE d4 = fma(c4, zRe, d);
            FP_TYPE e4 = fma(d4, zRe, e);
#endif
#if FP_PREC == 64
            FP_TYPE a4 = a;
            FP_TYPE b4 = a4*zRe + b;
            FP_TYPE c4 = b4*zRe + c;
            FP_TYPE d4 = c4*zRe + d;
            FP_TYPE e4 = d4*zRe + e;
#endif
            quarticSolve(a4, b4, c4, d4, e4, rootsRe, rootsIm);
        }
    }

    // Do one polish step, using Newton iteration,
    // but only do so if we do not screw up the roots.
    for (int i=0; i<5; i++) {
        FP_TYPE closestDist = closestRootDistance(i, rootsRe, rootsIm);
        if (rootsIm[i] == 0.0) {
            // We have a real root, use the more
            // efficient real Newton iteration.
            FP_TYPE fz = hornerReal(a, b, c, d, e, f, rootsRe[i]);
            FP_TYPE dfz = dhornerReal(da, db, dc, dd, de, rootsRe[i]);
            FP_TYPE corr = fz/dfz;
            if (isfinite(corr) && abs(corr) < 0.001*min(closestDist, abs(rootsRe[i]))) {
                rootsRe[i] -= corr;
            }
        }
        else {
            // We have a complex root, use the general
            // complex Newton iteration method.
            FP_TYPE numRe, numIm, denRe, denIm, corrRe, corrIm;
            horner(a, b, c, d, e, f, rootsRe[i], rootsIm[i], &numRe, &numIm);
            dhorner(da, db, dc, dd, de, rootsRe[i], rootsIm[i], &denRe, &denIm);
            cdiv(numRe, numIm, denRe, denIm, &corrRe, &corrIm);
            if (isfinite(corrRe) && isfinite(corrIm) && 
                HYPOT(corrRe, corrIm) < 0.001*min(closestDist, HYPOT(rootsRe[i], rootsIm[i]))) {
                rootsRe[i] -= corrRe;
                rootsIm[i] -= corrIm;
                // Also replace the next root, which is the complex conjugate.
                rootsRe[i+1] = rootsRe[i];
                rootsIm[i+1] = -rootsIm[i];
            }
            i++;  // Skip the next root, it is the complex conjugate.
        }
    }
}




// Computes the value of the polynomial, specified by the
// given coefs at the complex point z = zRe+i*zIm. The
// result is placed in the reference value result.
// Horner's multiplication scheme is used for efficiency.
static void horner(FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE d, FP_TYPE e, FP_TYPE f, FP_TYPE zRe, FP_TYPE zIm, FP_TYPE *resRe, FP_TYPE *resIm) {
    FP_TYPE r = 2*zRe;
    FP_TYPE s = zRe*zRe + zIm*zIm;
    FP_TYPE aa = a, bb = b;

#if FP_PREC == 53
    {
        FP_TYPE tmp = fma(r, aa, bb);  // bb + r*aa;
        bb = fma(-s, aa, c);           // c - s*aa;
        aa = tmp;
    }
    {
        FP_TYPE tmp = fma(r, aa, bb);  // bb + r*aa;
        bb = fma(-s, aa, d);           // d - s*aa;
        aa = tmp;
    }
    {
        FP_TYPE tmp = fma(r, aa, bb);  // bb + r*aa;
        bb = fma(-s, aa, e);           // e - s*aa;
        aa = tmp;
    }
    {
        FP_TYPE tmp = fma(r, aa, bb);  // bb + r*aa;
        bb = fma(-s, aa, f);           // f - s*aa;
        aa = tmp;
    }

    *resRe = fma(zRe, aa, bb);         // zRe*aa + bb;
    *resIm = zIm*aa;
#endif
#if FP_PREC == 64
    {
        FP_TYPE tmp = bb + r*aa;
        bb = c - s*aa;
        aa = tmp;
    }
    {
        FP_TYPE tmp = bb + r*aa;
        bb = d - s*aa;
        aa = tmp;
    }
    {
        FP_TYPE tmp = bb + r*aa;
        bb = e - s*aa;
        aa = tmp;
    }
    {
        FP_TYPE tmp = bb + r*aa;
        bb = f - s*aa;
        aa = tmp;
    }

    *resRe = zRe*aa + bb;
    *resIm = zIm*aa;
#endif
}




// Evaluates the quantic for a real given z.
static FP_TYPE hornerReal(FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE d, FP_TYPE e, FP_TYPE f, FP_TYPE z) {
#if FP_PREC == 53
    FP_TYPE fz = fma(a, z, b);
    fz = fma(fz, z, c);
    fz = fma(fz, z, d);
    fz = fma(fz, z, e);
    fz = fma(fz, z, f);
#endif
#if FP_PREC == 64
    FP_TYPE fz = a*z + b;
    fz = fz*z + c;
    fz = fz*z + d;
    fz = fz*z + e;
    fz = fz*z + f;
#endif
    return fz;
}




// Computes the derivative of the polynomial, specified by the
// given coefs at the complex point z = zRe + i*zIm. The
// result is placed in the reference value result.
// Horner's multiplication scheme is used for efficiency.
static void dhorner(FP_TYPE da, FP_TYPE db, FP_TYPE dc, FP_TYPE dd, FP_TYPE de, FP_TYPE zRe, FP_TYPE zIm, FP_TYPE *resRe, FP_TYPE *resIm) {
    FP_TYPE r = 2*zRe;
    FP_TYPE s = zRe*zRe + zIm*zIm;
    FP_TYPE aa = da, bb = db;

#if FP_PREC == 53
    {
        FP_TYPE tmp = fma(r, aa, bb);  // bb + r*aa;
        bb = fma(-s, aa, dc);         // dc - s*aa;
        aa = tmp;
    }
    {
        FP_TYPE tmp = fma(r, aa, bb);  // bb + r*aa;
        bb = fma(-s, aa, dd);         // dd - s*aa;
        aa = tmp;
    }
    {
        FP_TYPE tmp = fma(r, aa, bb);  // bb + r*aa;
        bb = fma(-s, aa, de);         // de - s*aa;
        aa = tmp;
    }

    *resRe = fma(zRe, aa, bb);     // zRe*aa + bb;
    *resIm = zIm*aa;
#endif
#if FP_PREC == 64
    {
        FP_TYPE tmp = bb + r*aa;
        bb = dc - s*aa;
        aa = tmp;
    }
    {
        FP_TYPE tmp = bb + r*aa;
        bb = dd - s*aa;
        aa = tmp;
    }
    {
        FP_TYPE tmp = bb + r*aa;
        bb = de - s*aa;
        aa = tmp;
    }

    *resRe = zRe*aa + bb;
    *resIm = zIm*aa;
#endif
}




// Computes the derivative of the quintic for a real given z.
static FP_TYPE dhornerReal(FP_TYPE da, FP_TYPE db, FP_TYPE dc, FP_TYPE dd, FP_TYPE de, FP_TYPE z) {
#if FP_PREC == 53
    FP_TYPE dfz = fma(da, z, db);
    dfz = fma(dfz, z, dc);
    dfz = fma(dfz, z, dd);
    dfz = fma(dfz, z, de);
#endif
#if FP_PREC == 64
    FP_TYPE dfz = da*z + db;
    dfz = dfz*z + dc;
    dfz = dfz*z + dd;
    dfz = dfz*z + de;
#endif
    return dfz;
}




// Complex division, using Smith's formula to avoid overflow or underflow.
static void cdiv(FP_TYPE numRe, FP_TYPE numIm, FP_TYPE denRe, FP_TYPE denIm, FP_TYPE *resRe, FP_TYPE *resIm) {
    FP_TYPE are = abs(denRe);
    FP_TYPE aim = abs(denIm);
    if (aim < are) {
        FP_TYPE doc = denIm / denRe;
        *resRe = (numRe + numIm * doc) / (denRe + denIm * doc);
        *resIm = (numIm - numRe * doc) / (denRe + denIm * doc);
    } 
    else {
        FP_TYPE cod = denRe / denIm;
        *resRe = (numIm + numRe * cod) / (denIm + denRe * cod);
        *resIm = (-numRe + numIm * cod) / (denIm + denRe * cod);
    }
}




// Use interval bisection to find a real root of the quintic equation.
// This always succeeds, because we can easily compute an interval,
// which contains an odd number of roots and then we can split the
// interval by two until it does not become smaller, due to limits
// in the system's precision. We cannot tell beforehand, which real
// solution is obtained, if more than one such solution exists.
static FP_TYPE bisect(FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE d, FP_TYPE e, FP_TYPE f) {
    if (a < 0) {
        a = -a; b=-b; c=-c; d=-d; e=-e; f=-f;
    }
    FP_TYPE maxRoot = abs(b/a);
    maxRoot = max(maxRoot, SQRT(abs(c/a)));
    maxRoot = max(maxRoot, CBRT(abs(d/a)));
    maxRoot = max(maxRoot, SQRT(SQRT(abs(e/a))));
    maxRoot = max(maxRoot, POW(abs(f/(2*a)), 0.2));
    FP_TYPE left = -maxRoot;
    while (hornerReal(a, b, c, d, e, f, left) > 0) {
        left += left;
    }
    FP_TYPE right = maxRoot;
    while (hornerReal(a, b, c, d, e, f, right) < 0) {
        right += right;
    }
    int niter = 0;
    while (niter < 500) {
        FP_TYPE mid = (left + right) * 0.5;
        if (left == mid || right == mid) {
           break;
        }
        FP_TYPE leftFz = hornerReal(a, b, c, d, e, f, mid);
        if (leftFz < 0) {
            left = mid;
        }
        else {
            right = mid;
        }
        niter++;
    }
    return (left + right) * 0.5;
}
    
    
    
    
static FP_TYPE closestRootDistance(int idx, FP_TYPE *rootsRe, FP_TYPE *rootsIm) {
    FP_TYPE distance = MAX_VALUE;
    for (int i=0; i<5; i++) {
        if (i == idx) {
            continue;
        }
        FP_TYPE dist = HYPOT(rootsRe[idx] - rootsRe[i], rootsIm[idx] - rootsIm[i]);
        if (dist < distance) {
            distance = dist;
        }
    }
    return distance;
}




// Like Java's Math.getExponent() method.
static inline int getExponent(FP_TYPE a) {
    int exp;
    frexp(a, &exp);
    return exp - 1;
}

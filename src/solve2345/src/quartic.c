#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "solve_utils.h"

#if FP_PREC == 53
#define SCALE1 0x1p200
#define SCALE2 0x1p400
#define SCALE3 0x1p600
#define SCALE4 0x1p800
#endif

#if FP_PREC == 64
#define SCALE1 0x1p1000L
#define SCALE2 0x1p2000L
#define SCALE3 0x1p3000L
#define SCALE4 0x1p4000L
#endif
/**
 *
 * @author woelen
 */
    

 
static void quadratic(FP_TYPE a, FP_TYPE b, FP_TYPE *zeror, FP_TYPE *zeroi, int rootOffset);
static void complexDecomposition_l1(FP_TYPE l1, FP_TYPE gamma, FP_TYPE l3, FP_TYPE gammal2, FP_TYPE *zeror, FP_TYPE *zeroi);
static void complexDecomposition(FP_TYPE gamma, FP_TYPE l3, FP_TYPE gammal2, FP_TYPE *zeror, FP_TYPE *zeroi);
static FP_TYPE cubic_B_shift_4args(FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE d);
static FP_TYPE cubic_B_shift_3args(FP_TYPE b, FP_TYPE c, FP_TYPE d);

#if FP_PREC == 64
static inline FP_TYPE hypotl_impl(FP_TYPE a, FP_TYPE b) {return SQRT(a*a+b*b);} 
#endif
    
/**
 * Solves a quartic equation with real coefficients of the form a*x^4 + b*x^3 + c*x^2 + d*x + e == 0.
 * The complete solution is returned in the supplied arrays, in no particular order.
 * @param a The coefficient of the quartic term
 * @param b The coefficient of the cubic term
 * @param c The coefficient of the quadratic term
 * @param d The coefficient of the linear term
 * @param e The constant term
 * @param zeror This is used for output. The first four elements of this array are
 * overwritten with the real part of the four zeros of the equation.
 * @param zeroi This is used for output. The first four elements of this array are
 * overwritten with the imaginary part of the four zeros of the equation.
 */
void quarticSolve(FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE d, FP_TYPE e, FP_TYPE *zeror, FP_TYPE *zeroi) {
    if (a == 1.0) {
        quarticSolveMonic(b, c, d, e, zeror, zeroi);
    }
    else {
        quarticSolveMonic(b/a, c/a, d/a, e/a, zeror, zeroi);
    }
}
    
    
    
    
    
//===================================================================
//---------------- The Low-Rank LDL^T Quartic Solver ----------------
//===================================================================
//          This program solves a quartic equation of the form:
//          Q(z) = z^4 + a*z^3 + b*z^2 + c*z + d 
//               = (z - z[0])*(z - z[1])*(z - z[2])*(z - z[3])
//          Input:  a,b,c,d  (The quartic coefficients)
//          Output: zeror, zeroi (arrays of 4 zeros)
//====================================================================

// Solves a monic quartic equation, using the LDL^T decomposition, as
// introduced by Peter Strobach. The original algorithm is modified
// to make it more robust and more accurate.
// The case d == 0 is handled by solving a cubic equation, it cannot
// be handled by the LDL^T method.


/**
 * Solves a quartic equation with real coefficients of the form x^4 + a*x^3 + b*x^2 + c*x + d == 0.
 * The complete solution is returned in the supplied arrays, in no particular order.
 * @param a The coefficient of the cubic term
 * @param b The coefficient of the quadratic term
 * @param c The coefficient of the linear term
 * @param d The constant term
 * @param zeror This is used for output. The first four elements of this array are
 * overwritten with the real part of the four zeros of the equation.
 * @param zeroi This is used for output. The first four elements of this array are
 * overwritten with the imaginary part of the four zeros of the equation.
 */
void quarticSolveMonic(FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE d, FP_TYPE *zeror, FP_TYPE *zeroi) {

    // Handle the case that the constant term equals 0 in a special way.
    // Solve a cubic equation, and add a root with the value 0.
    if (d == 0.0) {
        cubicSolveMonic(a, b, c, zeror, zeroi);
        zeror[3] = zeroi[3] = 0.0;
        return;
    }

    // Here starts the general LDL^T algorithm.


    // calculate the antidiagonal shift phi0.
    FP_TYPE phi0 = cubic_B_shift_4args(a, b, c, d);
    boolean scaled = false;
    boolean scaled2 = false;
    if (!isfinite(phi0)) {
        a = a/SCALE1;
        b = b/SCALE2;
        c = c/SCALE3;
        d = d/SCALE4;
        scaled = true;
        phi0 = cubic_B_shift_4args(a, b, c, d);
        if (!isfinite(phi0)) {
            a = a/SCALE1;
            b = b/SCALE2;
            c = c/SCALE3;
            d = d/SCALE4;
            scaled2 = true;
            phi0 = cubic_B_shift_4args(a, b, c, d);
        }
    }

    // Compute LDL^T parameters l1, l2, gamma, and gamma*l2,
    // together with their accuracies. Also determine, whether
    // a product of two real quadratic equations is determined,
    // or a product of two conjugate complex quadratic equations.
    FP_TYPE l1, l3, l3sumabs;
    FP_TYPE bmin2l3, bmin2l3sumabs;
    FP_TYPE gammal2, gamma;
    FP_TYPE gammal2Accuracy, gammaAccuracy;
    boolean realDecomposition;
    {
        l1 = a / 2;
        l3 = b / 6 + phi0 / 2;
        l3sumabs = abs(b/6) + abs(phi0/2);

        FP_TYPE al3 = a * l3;

        // Compute d2
        FP_TYPE l1l1 = l1*l1;
        FP_TYPE b23 = 2*b/3;
        bmin2l3 = b23 - phi0;
        bmin2l3sumabs = abs(b23) + abs(phi0);
        FP_TYPE d2 = bmin2l3 - l1l1;
        FP_TYPE d2Accuracy = abs(d2)/(MIN_NORMAL + bmin2l3sumabs + l1l1);
        if (d2Accuracy < EPS) {
            d2 = d2Accuracy = 0.0;
        }

        // Compute d2*l2
        FP_TYPE d2l2 = c - al3;
        FP_TYPE d2l2Accuracy = abs(d2l2)/(MIN_NORMAL + abs(c) + abs(a)*l3sumabs);
        d2l2 /= 2;

        // Compute d2*l2*l2
        FP_TYPE l3l3 = l3*l3;
        FP_TYPE d2l22 = d - l3l3;
        FP_TYPE d2l22Accuracy = abs(d2l22)/(MIN_NORMAL + abs(d) + l3sumabs*l3sumabs);
        if (d2l22Accuracy < EPS) {
            d2l22 = d2l22Accuracy = 0.0;
        }

        // The computation of the LDL^T parameters suffers from
        // numerical cancellation. Select the best possible option.
        if (d2 != 0.0 || d2l22 != 0.0) {
            if (d2Accuracy >= d2l22Accuracy) {
                // The value of d2 can be computed with the lowest 
                // numerical cancellation. Use that as the starting
                // point for computing the LDL^T parameters.
                realDecomposition = (d2 < 0);
                gamma = realDecomposition ? SQRT(-d2) : SQRT(d2);
                gammaAccuracy = d2Accuracy;
                FP_TYPE l2 = d2l2 / d2;      // Compute l2 from d2 and d2*l2.
                gammal2 = gamma * l2;
                gammal2Accuracy = min(d2l2Accuracy, d2Accuracy);
            }
            else {
                // We have the best-conditioned computation for d2*l2*l2.
                realDecomposition = (d2l22 < 0);
                gammal2 = realDecomposition
                             ? ((d2l2 < 0) ?  SQRT(-d2l22) : -SQRT(-d2l22))
                             : ((d2l2 < 0) ? -SQRT(d2l22)  :  SQRT(d2l22));
                gammal2Accuracy = d2l22Accuracy;
                gamma = abs(d2l2 / gammal2);
                gammaAccuracy = min(d2l2Accuracy, d2l22Accuracy);
            }
        }
        else {
            // In this situation, we have two FP_TYPE roots,
            // or we even have a single quadruple root.
            realDecomposition = true;
            gamma = gammal2 = 0.0;
            gammaAccuracy = gammal2Accuracy = 1.0;
        }
    }

    if (realDecomposition) {
        // We need a real decomposition.
        FP_TYPE aq = l1 + gamma, bq = l3 + gammal2;      // equations (23) and (24)
        FP_TYPE cq = l1 - gamma, dq = l3 - gammal2;      // equations (25) and (26)
        FP_TYPE bqdqAccuracy;   // Accuracy measure of both bq and dq.
        if (abs(dq) < abs(bq)) {
            dq = d / bq;  // Use of equation (30), case |b| > |d|
            bqdqAccuracy = abs(bq)/(MIN_NORMAL + l3sumabs + abs(gammal2)/gammal2Accuracy);
        } 
        else {
            bq = d / dq;  // Use of equation (30), case |b| <= |d|
            bqdqAccuracy = abs(dq)/(MIN_NORMAL + l3sumabs + abs(gammal2)/gammal2Accuracy);
        }

        FP_TYPE absaq = abs(aq);
        FP_TYPE abscq = abs(cq);
        FP_TYPE aqAccuracy = 0.0, cqAccuracy = 0.0;
        if (gamma > 0 && (absaq > 0 || abscq > 0)) {
            if (abscq < absaq) {
                aqAccuracy = absaq/(MIN_NORMAL + abs(l1) + gamma/gammaAccuracy);

                // We have three possible means of computing c,
                // as described in section 3.4, case 1. We take the
                // one with the least numerical cancellation.
                FP_TYPE c1 = a - aq;
                FP_TYPE c1Accuracy = abs(c1)/(abs(a) + absaq/aqAccuracy);
                if (c1Accuracy < 3*EPS) {
                    c1 = c1Accuracy = 0.0;
                }
                FP_TYPE c2Num = bmin2l3;
                FP_TYPE c2NumAccuracy = abs(c2Num)/(MIN_NORMAL + bmin2l3sumabs);
                FP_TYPE c2Accuracy = min(c2NumAccuracy, aqAccuracy);
                if (c2Accuracy < 3*EPS) {
                    c2Num = c2Accuracy = 0.0;
                }
                FP_TYPE aqdq = aq*dq;
                FP_TYPE c3Num = c - aqdq;
                FP_TYPE c3NumAccuracy = abs(c3Num)/(MIN_NORMAL + abs(c) + abs(aqdq)/min(aqAccuracy, bqdqAccuracy));
                FP_TYPE c3Accuracy = min(c3NumAccuracy, bqdqAccuracy);
                if (c3Accuracy < 3*EPS) {
                    c3Num = c3Accuracy = 0.0;
                }

                if (c1Accuracy >= c2Accuracy && c1Accuracy >= c3Accuracy) {
                    cq = c1;
                    cqAccuracy = c1Accuracy;
                }
                else if (c2Accuracy >= c1Accuracy && c2Accuracy >= c3Accuracy) {
                    cq = c2Num / aq;    // aq here is not 0 and it is accurate
                    cqAccuracy = c2Accuracy;
                }
                else {
                    cq = c3Num / bq;    // bq never is 0 and it is accurate
                    cqAccuracy = c3Accuracy;
                }
            }
            else {
                cqAccuracy = abscq/(MIN_NORMAL + abs(l1) + gamma/gammaAccuracy);

                // We have three possible means of computing a,
                // as described in section 3.4, case 2. We take the
                // one with the least numerical cancellation.
                FP_TYPE a1 = a - cq;
                FP_TYPE a1Accuracy = abs(a1)/(abs(a) + abscq/cqAccuracy);
                if (a1Accuracy < 3*EPS) {
                    a1 = a1Accuracy = 0.0;
                }
                FP_TYPE a2Num = bmin2l3;
                FP_TYPE a2NumAccuracy = abs(a2Num)/(MIN_NORMAL + bmin2l3sumabs);
                FP_TYPE a2Accuracy = min(a2NumAccuracy, cqAccuracy);
                if (a2Accuracy < 3*EPS) {
                    a2Num = a2Accuracy = 0.0;
                }
                FP_TYPE bqcq = bq*cq;
                FP_TYPE a3Num = c - bqcq;
                FP_TYPE a3NumAccuracy = abs(a3Num)/(MIN_NORMAL + abs(c) + abs(bqcq)/min(cqAccuracy, bqdqAccuracy));
                FP_TYPE a3Accuracy = min(a3NumAccuracy, bqdqAccuracy);
                if (a3Accuracy < 3*EPS) {
                    a3Num = a3Accuracy = 0.0;
                }

                if (a1Accuracy >= a2Accuracy && a1Accuracy >= a3Accuracy) {
                    aq = a1;
                    aqAccuracy = a1Accuracy;
                }
                else if (a2Accuracy >= a1Accuracy && a2Accuracy >= a3Accuracy) {
                    aq = a2Num / cq;    // cq here is not 0 and is accurate
                    aqAccuracy = a2Accuracy;
                }
                else {
                    aq = a3Num / dq;    // dq never is 0 and it is accurate
                    aqAccuracy = a3Accuracy;
                }
            }
        }

        // If the accuracy of (bq, dq) is not very good, then we use
        // an alternative way of computing these from the two equations,
        // mentioned below. The accuracy of (bq, dq) can be very low
        // if l3 is computed with large rounding error. This occurs,
        // when the original quartic has roots with a large dynamic
        // range. In such cases, (aq, cq) is computed with good accuracy,
        // while (bq, dq) has low accuracy. In such cases we derive
        // (bq, dq) from (aq, cq) and the original coefficients of the
        // quartic expression. This is described in section 3.5.
        if (bqdqAccuracy < min(aqAccuracy, cqAccuracy))
        {
            // eq. A: cq*bq^2 - c*bq + aq*d == 0  --> solve for bq
            // eq. B: aq*dq^2 - c*dq + cq*d == 0  --> solve for dq
            FP_TYPE aqcq = aq*cq;
            FP_TYPE aqcqd4 = 4*aqcq*d;
            FP_TYPE cc = c*c;
            FP_TYPE discr = cc - aqcqd4;
            if (discr > 0.0) {
                FP_TYPE sq = SQRT(discr);
                FP_TYPE discrAccuracy = discr/(MIN_NORMAL+cc+abs(aqcqd4)/min(aqAccuracy, cqAccuracy));
                FP_TYPE bq1, bq2, dq1, dq2;
                FP_TYPE csqAccuracy;
                if (c >= 0) {
                    FP_TYPE cplussq = c + sq;
                    bq1 = cplussq/(2*cq);     // use eq. A largest magnitude
                    dq1 = d/bq1;              // use D = b*d
                    dq2 = cplussq/(2*aq);     // use eq. B largest magnitude
                    bq2 = d/dq2;              // use D = b*d
                    csqAccuracy = cplussq/(MIN_NORMAL + c + sq/discrAccuracy);
                }
                else {
                    FP_TYPE cminussq = c - sq;
                    bq1 = cminussq/(2*cq);    // use eq. A largest magnitude
                    dq1 = d/bq1;              // use D = b*d
                    dq2 = cminussq/(2*aq);    // use eq. B largest magnitude
                    bq2 = d/dq2;              // use D = b*d
                    csqAccuracy = -cminussq/(MIN_NORMAL - c + sq/discrAccuracy);
                }
                if (csqAccuracy > bqdqAccuracy) {
                    FP_TYPE b_aqcq = b - aqcq;
                    FP_TYPE diff1 = b_aqcq - bq1 - dq1;
                    FP_TYPE diff2 = b_aqcq - bq2 - dq2;
                    if (abs(diff1) < abs(diff2)) {
                        bq = bq1;
                        dq = dq1;
                    }
                    else {
                        bq = bq2;
                        dq = dq2;
                    }
                }
            }
        }

        // roots of the quadratic factors:
        quadratic(aq, bq, zeror, zeroi, 0);
        quadratic(cq, dq, zeror, zeroi, 2);
    }
    else {
        // We need a complex decomposition.
        // We solve the following equations: 
        //          x^2 + (l1 + i*gamma)*x + l3 + i*gamma*l2 == 0
        //          x^2 + (l1 - i*gamma)*x + l3 - i*gamma*l2 == 0.
        complexDecomposition_l1(l1, gamma, l3, gammal2, zeror, zeroi);
    }

    // Undo any scaling, if we did that.
    if (scaled) {
        zeror[0] *= SCALE1;  zeroi[0] *= SCALE1;
        zeror[1] *= SCALE1;  zeroi[1] *= SCALE1;
        zeror[2] *= SCALE1;  zeroi[2] *= SCALE1;
        zeror[3] *= SCALE1;  zeroi[3] *= SCALE1;
        if (scaled2) {
            zeror[0] *= SCALE1;  zeroi[0] *= SCALE1;
            zeror[1] *= SCALE1;  zeroi[1] *= SCALE1;
            zeror[2] *= SCALE1;  zeroi[2] *= SCALE1;
            zeror[3] *= SCALE1;  zeroi[3] *= SCALE1;
        }
    }
}






/**
 * Solves a quartic equation with real coefficients of the form x^4 + b*x^2 + c*x + d == 0.
 * The complete solution is returned in the supplied arrays, in no particular order.
 * @param b The coefficient of the quadratic term
 * @param c The coefficient of the linear term
 * @param d The constant term
 * @param zeror This is used for output. The first four elements of this array are
 * overwritten with the real part of the four zeros of the equation.
 * @param zeroi This is used for output. The first four elements of this array are
 * overwritten with the imaginary part of the four zeros of the equation.
 */
void quarticSolveDepressed(FP_TYPE b, FP_TYPE c, FP_TYPE d, FP_TYPE *zeror, FP_TYPE *zeroi) {

    // Handle the case that the constant term equals 0 in a special way.
    // Solve a cubic equation, and add a root with the value 0.
    if (d == 0.0) {
        cubicSolveDepressed(b, c, zeror, zeroi);
        zeror[3] = zeroi[3] = 0.0;
        return;
    }

    // Here starts the general LDL^T algorithm.


    // calculate the antidiagonal shift phi0.
    FP_TYPE phi0 = cubic_B_shift_3args(b, c, d);
    boolean scaled = false;
    boolean scaled2 = false;
    if (!isfinite(phi0)) {
        b = b/SCALE2;
        c = c/SCALE3;
        d = d/SCALE4;
        scaled = true;
        phi0 = cubic_B_shift_3args(b, c, d);
        if (!isfinite(phi0)) {
            b = b/SCALE2;
            c = c/SCALE3;
            d = d/SCALE4;
            scaled2 = true;
            phi0 = cubic_B_shift_3args(b, c, d);
        }
    }

    // Compute LDL^T parameters l3, gamma, and gamma*l2,
    // and the accuracies of gamma and gamma*l2.
    FP_TYPE l3, gamma, gammal2;
    FP_TYPE gammaAccuracy, gammal2Accuracy;
    FP_TYPE l3sumabs;
    FP_TYPE bmin2l3, bmin2l3sumabs;
    boolean realDecomposition;
    {
        l3 = b / 6 + phi0 / 2;
        l3sumabs = abs(b/6) + abs(phi0/2);

        // Compute d2
        FP_TYPE b23 = 2*b/3;
        bmin2l3 = b23 - phi0;
        bmin2l3sumabs = abs(b23) + abs(phi0);
        FP_TYPE d2 = bmin2l3;
        FP_TYPE d2Accuracy = abs(d2)/(MIN_NORMAL + bmin2l3sumabs);
        if (d2Accuracy < EPS) {
            d2 = d2Accuracy = 0.0;
        }

        // Compute d2*l2
        FP_TYPE d2l2 = c/2;
        FP_TYPE d2l2Accuracy = (c != 0.0) ? 1.0 : 0.0;

        // Compute d2*l2*l2
        FP_TYPE l3l3 = l3*l3;
        FP_TYPE d2l22 = d - l3l3;
        FP_TYPE d2l22Accuracy = abs(d2l22)/(MIN_NORMAL + abs(d) + l3sumabs*l3sumabs);
        if (d2l22Accuracy < EPS) {
            d2l22 = d2l22Accuracy = 0.0;
        }

        // The computation of the LDL^T parameters suffers from
        // numerical cancellation. Select the best possible option.
        if (d2 != 0.0 || d2l22 != 0.0) {
            if (d2Accuracy >= d2l22Accuracy) {
                // The value of d2 can be computed with the lowest 
                // numerical cancellation. Use that as the starting
                // point for computing the LDL^T parameters.
                realDecomposition = (d2 < 0);
                gamma = realDecomposition ? SQRT(-d2) : SQRT(d2);
                gammaAccuracy = d2Accuracy;
                FP_TYPE l2 = d2l2 / d2;      // Compute l2 from d2 and d2*l2.
                gammal2 = gamma * l2;
                gammal2Accuracy = min(d2l2Accuracy, d2Accuracy);
            }
            else {
                // We have the best-conditioned computation for d2*l2*l2.
                realDecomposition = (d2l22 < 0);
                gammal2 = realDecomposition
                             ? ((d2l2 < 0) ?  SQRT(-d2l22) : -SQRT(-d2l22))
                             : ((d2l2 < 0) ? -SQRT(d2l22)  :  SQRT(d2l22));
                gammal2Accuracy = d2l22Accuracy;
                gamma = abs(d2l2 / gammal2);
                gammaAccuracy = min(d2l2Accuracy, d2l22Accuracy);
            }
        }
        else {
            // In this situation, we have two FP_TYPE roots,
            // or we even have a single quadruple root.
            realDecomposition = true;
            gamma = gammal2 = 0.0;
            gammaAccuracy = gammal2Accuracy = 1.0;
        }
    }

    if (realDecomposition) {
        // We need a real decomposition.
        // We solve the following equations: 
        //          x^2 + gamma*x + l3 + gamma*l2 == 0
        //          x^2 - gamma*x + l3 - gamma*l2 == 0.
        FP_TYPE aq =  gamma, bq = l3 + gammal2;      // equations (23) and (24)
        FP_TYPE cq = -gamma, dq = l3 - gammal2;      // equations (25) and (26)
        FP_TYPE bqdqAccuracy;   // Accuracy measure of both bq and dq.
        if (abs(dq) < abs(bq)) {
            dq = d / bq;  // Use of equation (30), case |b| > |d|
            bqdqAccuracy = abs(bq)/(MIN_NORMAL + l3sumabs + abs(gammal2)/gammal2Accuracy);
        } 
        else {
            bq = d / dq;  // Use of equation (30), case |b| <= |d|
            bqdqAccuracy = abs(dq)/(MIN_NORMAL + l3sumabs + abs(gammal2)/gammal2Accuracy);
        }

        FP_TYPE aqAccuracy = gammaAccuracy;

        // If the accuracy of (bq, dq) is not very good, then we use
        // an alternative way of computing these from the two equations,
        // mentioned below. The accuracy of (bq, dq) can be very low
        // if l3 is computed with large rounding error. This occurs,
        // when the original quartic has roots with a large dynamic
        // range. In such cases, (aq, cq) is computed with good accuracy,
        // while (bq, dq) has low accuracy. In such cases we derive
        // (bq, dq) from (aq, cq) and the original coefficients of the
        // quartic expression. This is described in section 3.5.
        if (bqdqAccuracy < aqAccuracy)
        {
            // eq. A: cq*bq^2 - c*bq + aq*d == 0  --> solve for bq
            // eq. B: aq*dq^2 - c*dq + cq*d == 0  --> solve for dq
            FP_TYPE aqcq = aq*cq;
            FP_TYPE aqcqd4 = 4*aqcq*d;
            FP_TYPE cc = c*c;
            FP_TYPE discr = cc - aqcqd4;
            if (discr > 0.0) {
                FP_TYPE sq = SQRT(discr);
                FP_TYPE discrAccuracy = discr/(MIN_NORMAL+cc+abs(aqcqd4)/aqAccuracy);
                FP_TYPE bq1, bq2, dq1, dq2;
                FP_TYPE csqAccuracy;
                if (c >= 0) {
                    FP_TYPE cplussq = c + sq;
                    bq1 = cplussq/(2*cq);     // use eq. A largest magnitude
                    dq1 = d/bq1;              // use D = b*d
                    dq2 = cplussq/(2*aq);     // use eq. B largest magnitude
                    bq2 = d/dq2;              // use D = b*d
                    csqAccuracy = cplussq/(MIN_NORMAL + c + sq/discrAccuracy);
                }
                else {
                    FP_TYPE cminussq = c - sq;
                    bq1 = cminussq/(2*cq);    // use eq. A largest magnitude
                    dq1 = d/bq1;              // use D = b*d
                    dq2 = cminussq/(2*aq);    // use eq. B largest magnitude
                    bq2 = d/dq2;              // use D = b*d
                    csqAccuracy = -cminussq/(MIN_NORMAL - c + sq/discrAccuracy);
                }
                if (csqAccuracy > bqdqAccuracy) {
                    FP_TYPE b_aqcq = b - aqcq;
                    FP_TYPE diff1 = b_aqcq - bq1 - dq1;
                    FP_TYPE diff2 = b_aqcq - bq2 - dq2;
                    if (abs(diff1) < abs(diff2)) {
                        bq = bq1;
                        dq = dq1;
                    }
                    else {
                        bq = bq2;
                        dq = dq2;
                    }
                }
            }
        }

        // roots of the quadratic factors:
        quadratic(aq, bq, zeror, zeroi, 0);
        quadratic(cq, dq, zeror, zeroi, 2);
    }
    else {
        // We need a complex decomposition.
        // We solve the following equations: 
        //          x^2 + i*gamma*x + l3 + i*gamma*l2 == 0
        //          x^2 - i*gamma*x + l3 - i*gamma*l2 == 0.
        complexDecomposition(gamma, l3, gammal2, zeror, zeroi);
    }

    // Undo any scaling, if we did that.
    if (scaled) {
        zeror[0] *= SCALE1;  zeroi[0] *= SCALE1;
        zeror[1] *= SCALE1;  zeroi[1] *= SCALE1;
        zeror[2] *= SCALE1;  zeroi[2] *= SCALE1;
        zeror[3] *= SCALE1;  zeroi[3] *= SCALE1;
        if (scaled2) {
            zeror[0] *= SCALE1;  zeroi[0] *= SCALE1;
            zeror[1] *= SCALE1;  zeroi[1] *= SCALE1;
            zeror[2] *= SCALE1;  zeroi[2] *= SCALE1;
            zeror[3] *= SCALE1;  zeroi[3] *= SCALE1;
        }
    }
}





/**
 * Solves a quartic equation with real coefficients of the form x^4 + c*x + d == 0.
 * The complete solution is returned in the supplied arrays, in no particular order.
 * @param c The coefficient of the linear term
 * @param d The constant term
 * @param zeror This is used for output. The first four elements of this array are
 * overwritten with the real part of the four zeros of the equation.
 * @param zeroi This is used for output. The first four elements of this array are
 * overwritten with the imaginary part of the four zeros of the equation.
 */
void quarticSolveFullyDepressed(FP_TYPE c, FP_TYPE d, FP_TYPE *zeror, FP_TYPE *zeroi) {

    // Handle the case that the constant term equals 0 in a special way.
    // Solve a cubic equation, and add a root with the value 0.
    if (d == 0.0) {
        FP_TYPE sqrt3 = SQRT(3.0);
        FP_TYPE root = CBRT(-c);
        zeror[0] = root;
        zeroi[0] = 0.0;
        zeror[2] = zeror[1] = -0.5*root;
        zeroi[2] = -(zeroi[1] = 0.5*sqrt3*root);
        zeror[3] = zeroi[3] = 0.0;
        return;
    }

    // Here starts the general LDL^T algorithm.

    // calculate the antidiagonal shift phi0.
    FP_TYPE phi0 = cubicSolveDepressedDominant(-4.0*d, -c*c);
    boolean scaled = false;
    boolean scaled2 = false;
    if (!isfinite(phi0)) {
        c = c/SCALE3;
        d = d/SCALE4;
        scaled = true;
        phi0 = cubicSolveDepressedDominant(-4.0*d, -c*c);
        if (!isfinite(phi0)) {
            c = c/SCALE3;
            d = d/SCALE4;
            scaled2 = true;
            phi0 = cubicSolveDepressedDominant(-4.0*d, -c*c);
        }
    }

    // Compute LDL^T parameters l3, gamma, and gamma*l2.
    // Also determine whether we need to solve a pair
    // of complex equations, or a pair of real equations.
    FP_TYPE l3, gamma, gammal2;
    boolean realDecomposition;
    {
        l3 = phi0 / 2;

        // Compute d2
        FP_TYPE d2 = -phi0;

        // Compute d2*l2
        FP_TYPE d2l2 = c/2;

        // Compute d2*l2*l2
        FP_TYPE l3l3 = l3*l3;
        FP_TYPE d2l22 = d - l3l3;

        // Compute the LDL^T parameters gamma and gamma*l2.
        if (d2 != 0.0 || d2l22 != 0.0) {
            if (l3l3 > 0.5*abs(d)) {
                realDecomposition = (d2 < 0);
                gamma = realDecomposition ? SQRT(-d2) : SQRT(d2);
                FP_TYPE l2 = d2l2 / d2;      // Compute l2 from d2 and d2*l2.
                gammal2 = gamma * l2;
            }
            else {
                // We have the best-conditioned computation for d2*l2*l2.
                realDecomposition = (d2l22 < 0);
                gammal2 = realDecomposition
                             ? ((d2l2 < 0) ?  SQRT(-d2l22) : -SQRT(-d2l22))
                             : ((d2l2 < 0) ? -SQRT(d2l22)  :  SQRT(d2l22));
                gamma = abs(d2l2 / gammal2);
            }
        }
        else {
            // In this situation, we have two pairs of real FP_TYPE roots.
            realDecomposition = true;
            gamma = gammal2 = 0.0;
        }
    }

    if (realDecomposition) {
        // We need a real decomposition.
        // We solve the following equations: 
        //          x^2 + gamma*x + l3 + gamma*l2 == 0
        //          x^2 - gamma*x + l3 - gamma*l2 == 0.
        FP_TYPE aq =  gamma, bq = l3 + gammal2;      // equations (23) and (24)
        FP_TYPE cq = -gamma, dq = l3 - gammal2;      // equations (25) and (26)
        if (abs(dq) < abs(bq)) {
            dq = d / bq;  // Use of equation (30), case |b| > |d|
        } 
        else {
            bq = d / dq;  // Use of equation (30), case |b| <= |d
        }

        // roots of the quadratic factors:
        quadratic(aq, bq, zeror, zeroi, 0);
        quadratic(cq, dq, zeror, zeroi, 2);
    }
    else {
        // We need a complex decomposition.
        // We solve the following equations: 
        //          x^2 + i*gamma*x + l3 + i*gamma*l2 == 0
        //          x^2 - i*gamma*x + l3 - i*gamma*l2 == 0.
        complexDecomposition(gamma, l3, gammal2, zeror, zeroi);
    }

    // Undo any scaling, if we did that.
    if (scaled) {
        zeror[0] *= SCALE1;  zeroi[0] *= SCALE1;
        zeror[1] *= SCALE1;  zeroi[1] *= SCALE1;
        zeror[2] *= SCALE1;  zeroi[2] *= SCALE1;
        zeror[3] *= SCALE1;  zeroi[3] *= SCALE1;
        if (scaled2) {
            zeror[0] *= SCALE1;  zeroi[0] *= SCALE1;
            zeror[1] *= SCALE1;  zeroi[1] *= SCALE1;
            zeror[2] *= SCALE1;  zeroi[2] *= SCALE1;
            zeror[3] *= SCALE1;  zeroi[3] *= SCALE1;
        }
    }
}








/**
 * Solves a quartic equation with real coefficients of the form a*x^4 + c*x^2 + e == 0.
 * The complete solution is returned in the supplied arrays, in no particular order.
 * @param a The coefficient of the quartic term
 * @param c The coefficient of the quadratic term
 * @param e The constant term
 * @param zeror This is used for output. The first four elements of this array are
 * overwritten with the real part of the four zeros of the equation.
 * @param zeroi This is used for output. The first four elements of this array are
 * overwritten with the imaginary part of the four zeros of the equation.
 */
void quarticSolveBiquadratic(FP_TYPE a, FP_TYPE c, FP_TYPE e, FP_TYPE *zeror, FP_TYPE *zeroi) {
    if (a == 1.0) {
        quadraticSolveMonic(c, e, zeror, zeroi);
    }
    else {
        quadraticSolve(a, c, e, zeror, zeroi);
    }

    if (zeroi[0] == 0.0) {
        // We have two real solutions of the quadratic equation. 
        // These solutions can each lead to imaginary roots or 
        // two real roots of opposite sign, depending on sign.
        FP_TYPE root1 = zeror[0];
        FP_TYPE root2 = zeror[1];
        if (root1 >= 0) {
            zeror[1] = SQRT(root1);
            zeror[0] = -zeror[1];
            // No need to set imaginary part of the roots,
            // this was done by the quadratic solver already.
        }
        else {
            zeror[0] = zeror[1] = 0.0;
            zeroi[0] = SQRT(-root1);
            zeroi[1] = -zeroi[0];
        }

        if (root2 >= 0) {
            zeror[3] = SQRT(root2);
            zeror[2] = -zeror[3];
            zeroi[2] = zeroi[3] = 0.0;
        }
        else {
            zeror[2] = zeror[3] = 0.0;
            zeroi[2] = SQRT(-root2);
            zeroi[3] = -zeroi[2];
        }
    }
    else {
        // We have two complex solutions of the quadratic equations.
        // Take one complex solution of the quadratic.
        FP_TYPE re = zeror[0];   // This can have any sign.
        FP_TYPE im = zeroi[0];   // quadraticSolve() assures this is positive.

        // Take the square root of the solution and
        // finally take all 4 possible combinations
        // of +/- sqRe and +/- sqIm. No need to 
        // guard for zero, we always have non-zero
        // values. Also no need to check the sign of
        // the imaginary part (that is necessary for
        // true root extraction), because we create
        // all possible combinations of +/- re and
        // +/- im.
        FP_TYPE r = HYPOT(re, im);
        FP_TYPE zre, zim;
        if (re >= 0) {
            FP_TYPE sq = SQRT(2*(r + re));
            zre = 0.5 * sq;
            zim = im / sq;
        } 
        else {
            FP_TYPE sq = SQRT(2*(r - re));
            zre = im / sq;
            zim = 0.5 * sq;
        }
        zeror[0] = -zre;  zeroi[0] =  zim;
        zeror[1] = -zre;  zeroi[1] = -zim;
        zeror[2] =  zre;  zeroi[2] =  zim;
        zeror[3] =  zre;  zeroi[3] = -zim;
    }
}



    
    
    
    
    
    
/**
 * Given a set of roots and the coefficient of the leading term, compute
 * all coefficients of the quartic equation.
 * @param factor The coefficient of the leading term (the quartic term).
 * @param nReal Tells whether the equation has 0, 2, or 4 real roots. The 
 * next four parameters have meaning, depending on this parameter. The
 * first nReal parameters are regarded as real roots, the rest of these
 * four parameters are regarded as the real and imaginary part of complex
 * conjugate roots.
 * @param a If nReal equals 2 or 4, then this is interpreted as one of the
 * real roots of the equation. If nReal equals 0, then this is interpreted
 * as the real part of a pair of complex conjugate roots.
 * @param b If nReal equals 2 or 4, then this is interpreted as one of the
 * real roots of the equation. If nReal equals 0, then this is interpreted
 * as the imaginary part of a pair of complex conjugate roots, the other
 * conjugate complex root has imaginary part equal to -b.
 * @param c If nReal equals 4, then this is interpreted as one of the real
 * roots of the equation. If nReal equals 0 or 2, then this is interpreted
 * as the real part of a pair of complex conjugate roots.
 * @param d If nReal equals 4, then this is interpreted as one of the real
 * roots of the equation. If nReal equals 0 or 2, then this is interpreted
 * as the imaginary part of a pair of complex conjugate roots, the other
 * conjugate complex root has imaginary part equal to -d.
 * @param coefs An output array, which is filled with the coefficients
 * of the quartic. This array must have a length at least equal to 5. 
 * At index 0, the constant term is stored, at index 4 the coefficient
 * of the quartic term is stored.
 */
void quarticCoefs(FP_TYPE factor, int nReal, FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE d, FP_TYPE *coefs) {
    // First, compute the quadratic factor x^2 + p*x + q, based on a and b.
    FP_TYPE p, q;
    switch (nReal) {
        case 0:
            // There are no real solutions, the arguments a and b are regarded
            // as real and imaginary parts of the complex roots a ± i*b.
            p = -2.0 * a;
            q = a*a + b*b;
            break;
        case 2:
        case 4:
            // Either 2, or 4 real roots. The parameters a and b are
            // regarded as real roots.
            p = -(a + b);
            q = a * b;
            break;
        default:
            fprintf(stderr, "Real quartic equation has either four real solutions, two real solutions or no real solution.\n");
            exit(1);
    }

    // Next, compute the quadratic factor, based on c and d.
    FP_TYPE r, s;
    switch (nReal) {
        case 0:
        case 2:
            // There are no real solutions or only two real solutions.
            // The arguments c and d are regarded as real and imaginary
            // parts of the complex roots c ± i*d.
            r = -2.0 * c;
            s = c*c + d*d;
            break;
        case 4:
            // There are 4 real roots. The parameters c and d are
            // regarded as real roots.
            r = -(c + d);
            s = c * d;
            break;
        default:
            fprintf(stderr, "Real quartic equation has either four real solutions, two real solutions or no real solution.\n");
            exit(1);
    }

    // Compute factor * (x^2 + px + q) * (x^2 + rx + s).
    coefs[0] = factor * q * s;
    coefs[1] = factor * (p*s + q*r);
    coefs[2] = factor * (p*r + q + s);
    coefs[3] = factor * (p + r);
    coefs[4] = factor;
}






// solve a quadratic equation 
static void quadratic(FP_TYPE a, FP_TYPE b, FP_TYPE *zeror, FP_TYPE *zeroi, int rootOffset) {
    FP_TYPE diskr = a * a - 4 * b;
    if (diskr >= 0.0) {
        FP_TYPE div, zmax, zmin;
        if (a >= 0.0) {
            div = -a - SQRT(diskr);
        } 
        else {
            div = -a + SQRT(diskr);
        }

        zmax = div / 2;
        if (zmax == 0.0) {
            zmin = 0.0;
        } 
        else {
            zmin = b / zmax;
        }

        zeror[rootOffset] = zmax;
        zeror[rootOffset+1] = zmin;
        zeroi[rootOffset] = zeroi[rootOffset+1] = 0.0;
    } 
    else {
        FP_TYPE sq = SQRT(-diskr)/2;
        zeror[rootOffset] = zeror[rootOffset+1] = -a/2;
        zeroi[rootOffset] = sq;
        zeroi[rootOffset+1] = -sq;
    }
}





static void complexDecomposition_l1(FP_TYPE l1, FP_TYPE gamma, FP_TYPE l3, FP_TYPE gammal2, FP_TYPE *zeror, FP_TYPE *zeroi) {
    FP_TYPE acxRe = l1, acxIm = gamma;          // equation (10)
    FP_TYPE bcxRe = l3, bcxIm = gammal2;        // equation (10)

    // Computation of complex discriminant acx*acx - 4*bcx
    // and the square root of that discriminant.
    FP_TYPE discrRe = (acxRe*acxRe - acxIm*acxIm) - 4*bcxRe;
    FP_TYPE discrIm = 2*acxRe*acxIm - 4*bcxIm;
    FP_TYPE r = HYPOT(discrRe, discrIm);
    FP_TYPE sqRe, sqIm;
    if (discrRe != 0.0 || discrIm != 0.0) {
        if (discrRe >= 0) {
            FP_TYPE s = SQRT(2*(r + discrRe)); 
            sqRe = 0.5 * s;
            sqIm = abs(discrIm) / s;
        } 
        else {
            FP_TYPE s = SQRT(2*(r - discrRe));
            sqRe = abs(discrIm) / s;
            sqIm = 0.5 * s;
        }
        if (discrIm < 0) {
            sqIm = -sqIm;
        }
    }
    else {
        sqRe = sqIm = 0.0;
    }


    FP_TYPE zx1Re = (-acxRe + sqRe)/2, zx1Im = (-acxIm + sqIm)/2;
    FP_TYPE zx2Re = (-acxRe - sqRe)/2, zx2Im = (-acxIm - sqIm)/2;
    FP_TYPE zx1RR = zx1Re*zx1Re + zx1Im*zx1Im;
    FP_TYPE zx2RR = zx2Re*zx2Re + zx2Im*zx2Im;

    // zxmin = bcx.div(zxmax)
    FP_TYPE zxMaxRe, zxMaxIm, zxMinRe, zxMinIm;
    if (zx1RR > zx2RR) {
        // zx1 has highest magnitude, and hence low cancellation.
        // Compute zxmin, with zx1 as zxmax.
        zxMaxRe = zx1Re;
        zxMaxIm = zx1Im;
        zxMinRe = (bcxRe*zx1Re + bcxIm*zx1Im) / zx1RR;
        zxMinIm = (bcxIm*zx1Re - bcxRe*zx1Im) / zx1RR;
    }
    else {
        // zx2 has highest magnitude, and hence low cancellation.
        // Compute zxmin, with zx2 as zxmax.
        zxMaxRe = zx2Re;
        zxMaxIm = zx2Im;
        zxMinRe = (bcxRe*zx2Re + bcxIm*zx2Im) / zx2RR;
        zxMinIm = (bcxIm*zx2Re - bcxRe*zx2Im) / zx2RR;
    }

    zeror[1] = zeror[0] = zxMinRe;
    zeroi[1] = -(zeroi[0] = zxMinIm);
    zeror[3] = zeror[2] = zxMaxRe;
    zeroi[3] = -(zeroi[2] = zxMaxIm);
}




static void complexDecomposition(FP_TYPE gamma, FP_TYPE l3, FP_TYPE gammal2, FP_TYPE *zeror, FP_TYPE *zeroi) {
    FP_TYPE acxIm = gamma;                      // equation (10)
    FP_TYPE bcxRe = l3, bcxIm = gammal2;        // equation (10)

    // Computation of complex discriminant acx*acx - 4*bcx
    // and the square root of that discriminant.
    FP_TYPE discrRe = -acxIm*acxIm - 4*bcxRe;
    FP_TYPE discrIm = -4*bcxIm;
    FP_TYPE r = HYPOT(discrRe, discrIm);
    FP_TYPE sqRe, sqIm;
    if (discrRe != 0.0 || discrIm != 0.0) {
        if (discrRe >= 0) {
            FP_TYPE s = SQRT(2*(r + discrRe)); 
            sqRe = 0.5 * s;
            sqIm = abs(discrIm) / s;
        } 
        else {
            FP_TYPE s = SQRT(2*(r - discrRe));
            sqRe = abs(discrIm) / s;
            sqIm = 0.5 * s;
        }
        if (discrIm < 0) {
            sqIm = -sqIm;
        }
    }
    else {
        sqRe = sqIm = 0.0;
    }


    FP_TYPE zx1Re =  sqRe/2, zx1Im = (-acxIm + sqIm)/2;
    FP_TYPE zx2Re = -sqRe/2, zx2Im = (-acxIm - sqIm)/2;

    // zxmin = bcx.div(zxmax)
    FP_TYPE zxMaxRe, zxMaxIm, zxMinRe, zxMinIm;
    if (abs(zx1Im) > abs(zx2Im)) {
        // zx1 has highest magnitude, and hence low cancellation.
        // Compute zxmin, with zx1 as zxmax.
        FP_TYPE zx1RR = zx1Re*zx1Re + zx1Im*zx1Im;
        zxMaxRe = zx1Re;
        zxMaxIm = zx1Im;
        zxMinRe = (bcxRe*zx1Re + bcxIm*zx1Im) / zx1RR;
        zxMinIm = (bcxIm*zx1Re - bcxRe*zx1Im) / zx1RR;
    }
    else {
        // zx2 has highest magnitude, and hence low cancellation.
        // Compute zxmin, with zx2 as zxmax.
        FP_TYPE zx2RR = zx2Re*zx2Re + zx2Im*zx2Im;
        zxMaxRe = zx2Re;
        zxMaxIm = zx2Im;
        zxMinRe = (bcxRe*zx2Re + bcxIm*zx2Im) / zx2RR;
        zxMinIm = (bcxIm*zx2Re - bcxRe*zx2Im) / zx2RR;
    }

    zeror[1] = zeror[0] = zxMinRe;
    zeroi[1] = -(zeroi[0] = zxMinIm);
    zeror[3] = zeror[2] = zxMaxRe;
    zeroi[3] = -(zeroi[2] = zxMaxIm);
}





// phi0 calculated as the dominant root of the B-shifted characteristic 
// depressed cubic using the parabola/reciprocal intersection method as
// described in Strobachs paper. 
// All in standard FP_TYPE precison. No extended precision arithmetic involved.
static FP_TYPE cubic_B_shift_4args(FP_TYPE a, FP_TYPE b, FP_TYPE c, FP_TYPE d) {

    FP_TYPE phi0, g, h, gg, hh, aq, bq, cq, dq, s, diskr;

    // the B-shift:
    diskr = 9 * a * a - 24 * b;  // equation (3.58) of Strobachs paper
    if (diskr > 0.0) {
        diskr = SQRT(diskr);
        if (a > 0.0) {
            s = -2 * b / (3 * a + diskr);   // equation (3.58) of Strobachs paper
        } 
        else {
            s = -2 * b / (3 * a - diskr);   // equation (3.58) of Strobachs paper
        }
    } 
    else {
        s = -a / 4;      // equations (3.59-60)
    }

    // the shift transformation (Horner forms):
    aq = a + 4 * s;                               // equation (3.45) of Strobachs paper
    bq = b + 3 * s * (a + 2 * s);                 // equation (3.46) of Strobachs paper
    cq = c + s * (2 * b + s * (3 * a + 4 * s));   // equation (3.47) of Strobachs paper
    dq = d + s * (c + s * (b + s * (a + s)));     // equation (3.48) of Strobachs paper

    gg = bq*bq/9;
    hh = aq*cq;
    g = hh - 4*dq - 3*gg;                             // equation (3.49) of Strobachs paper
    h = (8*dq + hh - 2*gg)*bq/3 - cq*cq - dq*aq*aq;   // equation (3.50) of Strobachs paper

    // call the parabola/reciprocal intersection method:      
    phi0 = cubicSolveDepressedDominant(g, h);
    return phi0;
}



// phi0 calculated as the dominant root of the B-shifted characteristic 
// depressed cubic using the parabola/reciprocal intersection method as
// described in Strobachs paper. 
// All in standard FP_TYPE precison. No extended precision arithmetic involved. 
static FP_TYPE cubic_B_shift_3args(FP_TYPE b, FP_TYPE c, FP_TYPE d) {

    if (b >= 0.0) {
        // We can use a strongly simplified mechanism for finding g and h,
        // because we simply need to take s equal to 0 in this case.
        FP_TYPE bb = b*b;
        return cubicSolveDepressedDominant(-4*d - bb/3, (8*d-2*bb/9)*b/3 - c*c);
    }

    FP_TYPE phi0, g, h, gg, hh, aq, bq, cq, dq, s;

    s = -SQRT(-b/6);

    // the shift transformation (Horner forms):
    aq = 4 * s;                               // equation (3.45) of Strobachs paper
    bq = b + 6 * s * s;                       // equation (3.46) of Strobachs paper
    cq = c + s * (2 * b + 4 * s * s);         // equation (3.47) of Strobachs paper
    dq = d + s * (c + s * (b + s * s));       // equation (3.48) of Strobachs paper

    gg = bq*bq/9;
    hh = aq*cq;
    g = hh - 4*dq - 3*gg;                             // equation (3.49) of Strobachs paper
    h = (8*dq + hh - 2*gg)*bq/3 - cq*cq - dq*aq*aq;   // equation (3.50) of Strobachs paper

    // call the parabola/reciprocal intersection method:      
    phi0 = cubicSolveDepressedDominant(g, h);
    return phi0;
}

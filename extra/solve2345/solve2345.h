#ifndef _SOLVE2345_H
#define _SOLVE2345_H


void quadraticSolve(double a, double b, double c, double *zeror, double *zeroi);
void quadraticSolveMonic(double a, double b, double *zeror, double *zeroi);
void quadraticCoefs(double factor, int nReal, double a, double b, double *coefs);

void cubicSolve(double a, double b, double c, double d, double *zeror, double *zeroi);
void cubicSolveMonic(double a, double b, double c, double *zeror, double *zeroi);
void cubicSolveDepressed(double g, double h, double *zeror, double *zeroi);
double cubicSolveDepressedDominant(double g, double h);
void cubicCoefs(double factor, int nReal, double a, double b, double c, double *coefs);

void quarticSolve(double a, double b, double c, double d, double e, double *zeror, double *zeroi);
void quarticSolveMonic(double a, double b, double c, double d, double *zeror, double *zeroi);
void quarticSolveDepressed(double b, double c, double d, double *zeror, double *zeroi);
void quarticSolveFullyDepressed(double c, double d, double *zeror, double *zeroi);
void quarticSolveBiquadratic(double a, double c, double e, double *zeror, double *zeroi);
void quarticCoefs(double factor, int nReal, double a, double b, double c, double d, double *coefs);

void quinticSolve(double a, double b, double c, double d, double e, double f,
                  double *zeror, double *zeroi);
void quinticSolveMonic(double a, double b, double c, double d, double e,
                       double *zeror, double *zeroi);
double quinticSolveBringJerrardDominant(double g, double h);
void quinticSolveBringJerrard(double g, double h, double *zeror, double *zeroi);
void quinticCoefs(double factor, int nReal, double a, double b, double c, double d, double e, double *coefs);




void quadraticSolvel(long double a, long double b, long double c, long double *zeror, long double *zeroi);
void quadraticSolveMonicl(long double a, long double b, long double *zeror, long double *zeroi);
void quadraticCoefsl(long double factor, int nReal, long double a, long double b, long double *coefs);

void cubicSolvel(long double a, long double b, long double c, long double d, long double *zeror, long double *zeroi);
void cubicSolveMonicl(long double a, long double b, long double c, long double *zeror, long double *zeroi);
void cubicSolveDepressedl(long double g, long double h, long double *zeror, long double *zeroi);
long double cubicSolveDepressedDominantl(long double g, long double h);
void cubicCoefsl(long double factor, int nReal, long double a, long double b, long double c, long double *coefs);

void quarticSolvel(long double a, long double b, long double c, long double d, long double e, long double *zeror, long double *zeroi);
void quarticSolveMonicl(long double a, long double b, long double c, long double d, long double *zeror, long double *zeroi);
void quarticSolveDepressedl(long double b, long double c, long double d, long double *zeror, long double *zeroi);
void quarticSolveFullyDepressedl(long double c, long double d, long double *zeror, long double *zeroi);
void quarticSolveBiquadraticl(long double a, long double c, long double e, long double *zeror, long double *zeroi);
void quarticCoefsl(long double factor, int nReal, long double a, long double b, long double c, long double d, long double *coefs);

void quinticSolvel(long double a, long double b, long double c, long double d, long double e, long double f,
                   long double *zeror, long double *zeroi);
void quinticSolveMonicl(long double a, long double b, long double c, long double d, long double e,
                        long double *zeror, long double *zeroi);
long double quinticSolveBringJerrardDominantl(long double g, long double h);
void quinticSolveBringJerrardl(long double g, long double h, long double *zeror, long double *zeroi);
void quinticCoefsl(long double factor, int nReal, long double a, long double b, long double c, long double d, long double e, long double *coefs);


#endif

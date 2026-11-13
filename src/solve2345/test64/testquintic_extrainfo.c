#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>

#include "solve2345.h"
#include "test64.h"




static void printinfo(char *fmt, long double err, long double a, long double b, long double c, long double d, long double e, long double *coefs) {
    printf(fmt, err);
    printf("   %.18Le  %.18Le  %.18Le  %.18Le  %.18Le\n", a, b, c, d, e);
    printf("   coefs = %.18Le  %.18Le  %.18Le  %.18Le  %.18Le  %.18Le\n", coefs[0], coefs[1], coefs[2], coefs[3], coefs[4], coefs[5]);
    printf("----------------\n");
}

/**
 *
 * @author woelen
 */
    
int main(int argc, char **argv) {
    srand48(time(NULL));
    long double coefs[6];
    long double zeror[5];
    long double zeroi[5];
    for (int i=0; i<1000000000; i++) {
        long double a = 5 - 10*drand48();
        long double b = 5 - 10*drand48();
        long double c = 5 - 10*drand48();
        long double d = 5 - 10*drand48();
        long double e = 5 - 10*drand48();
        int nReal = 1 + (lrand48() % 3) * 2;
        quinticCoefsl(0.1+drand48(), nReal, a, b, c, d, e, coefs);
        quinticSolvel(coefs[5], coefs[4], coefs[3], coefs[2], coefs[1], coefs[0], zeror, zeroi);
        switch (nReal) {
            case 1:
                if (relErrorReal(5, a, zeror, zeroi) > 1e-4) printinfo("1R %Le\n", relErrorReal(5, a, zeror, zeroi), a, b, c, d, e, coefs);
                if (relErrorComplex(5, b, c, zeror, zeroi) > 1e-4) printinfo("1C %Le\n", relErrorComplex(5, b, c, zeror, zeroi), a, b, c, d, e, coefs);
                if (relErrorComplex(5, b, -c, zeror, zeroi) > 1e-4) printinfo("1C %Le\n", relErrorComplex(5, b, -c, zeror, zeroi), a, b, c, d, e, coefs);
                if (relErrorComplex(5, d, e, zeror, zeroi) > 1e-4) printinfo("1C %Le\n", relErrorComplex(5, d, e, zeror, zeroi), a, b, c, d, e, coefs);
                if (relErrorComplex(5, d, -e, zeror, zeroi) > 1e-4) printinfo("1C %Le\n", relErrorComplex(5, d, -e, zeror, zeroi), a, b, c, d, e, coefs);
                break;
            case 3:
                if (relErrorReal(5, a, zeror, zeroi) > 1e-4) printinfo("3R %Le\n", relErrorReal(5, a, zeror, zeroi), a, b, c, d, e, coefs);
                if (relErrorReal(5, b, zeror, zeroi) > 1e-4) printinfo("3R %Le\n", relErrorReal(5, b, zeror, zeroi), a, b, c, d, e, coefs);
                if (relErrorReal(5, c, zeror, zeroi) > 1e-4) printinfo("3R %Le\n", relErrorReal(5, c, zeror, zeroi), a, b, c, d, e, coefs);
                if (relErrorComplex(5, d, e, zeror, zeroi) > 1e-4) printinfo("3C %Le\n", relErrorComplex(5, d, e, zeror, zeroi), a, b, c, d, e, coefs);
                if (relErrorComplex(5, d, -e, zeror, zeroi) > 1e-4) printinfo("3C %Le\n", relErrorComplex(5, d, -e, zeror, zeroi), a, b, c, d, e, coefs);
                break;
            case 5:
                if (relErrorReal(5, a, zeror, zeroi) > 1e-4) printinfo("5R %Le\n", relErrorReal(5, a, zeror, zeroi), a, b, c, d, e, coefs);
                if (relErrorReal(5, b, zeror, zeroi) > 1e-4) printinfo("5R %Le\n", relErrorReal(5, b, zeror, zeroi), a, b, c, d, e, coefs);
                if (relErrorReal(5, c, zeror, zeroi) > 1e-4) printinfo("5R %Le\n", relErrorReal(5, c, zeror, zeroi), a, b, c, d, e, coefs);
                if (relErrorReal(5, d, zeror, zeroi) > 1e-4) printinfo("5R %Le\n", relErrorReal(5, d, zeror, zeroi), a, b, c, d, e, coefs);
                if (relErrorReal(5, e, zeror, zeroi) > 1e-4) printinfo("5R %Le\n", relErrorReal(5, e, zeror, zeroi), a, b, c, d, e, coefs);
                break;
        }
    }
}

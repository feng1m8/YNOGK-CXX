#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>

#include "solve2345.h"
#include "test53.h"


/**
 *
 * @author woelen
 */
    
int main(int argc, char **argv) {
    srand48(time(NULL));
    double coefs[6];
    double zeror[5];
    double zeroi[5];
    for (int i=0; i<100000000; i++) {
        double a = 5 - 10*drand48();
        double b = 5 - 10*drand48();
        double c = 5 - 10*drand48();
        double d = 5 - 10*drand48();
        double e = 5 - 10*drand48();
        int nReal = 1 + (lrand48() % 3) * 2;
        quinticCoefs(0.1+drand48(), nReal, a, b, c, d, e, coefs);
        quinticSolve(coefs[5], coefs[4], coefs[3], coefs[2], coefs[1], coefs[0], zeror, zeroi);
        switch (nReal) {
            case 1:
                if (relErrorReal(5, a, zeror, zeroi) > 1e-5) printf("1R %le\n", relErrorReal(5, a, zeror, zeroi));
                if (relErrorComplex(5, b, c, zeror, zeroi) > 1e-5) printf("1C %le\n", relErrorComplex(5, b, c, zeror, zeroi));
                if (relErrorComplex(5, b, -c, zeror, zeroi) > 1e-5) printf("1C %le\n", relErrorComplex(5, b, -c, zeror, zeroi));
                if (relErrorComplex(5, d, e, zeror, zeroi) > 1e-5) printf("1C %le\n", relErrorComplex(5, d, e, zeror, zeroi));
                if (relErrorComplex(5, d, -e, zeror, zeroi) > 1e-5) printf("1C %le\n", relErrorComplex(5, d, -e, zeror, zeroi));
                break;
            case 3:
                if (relErrorReal(5, a, zeror, zeroi) > 1e-5) printf("3R %le\n", relErrorReal(5, a, zeror, zeroi));
                if (relErrorReal(5, b, zeror, zeroi) > 1e-5) printf("3R %le\n", relErrorReal(5, b, zeror, zeroi));
                if (relErrorReal(5, c, zeror, zeroi) > 1e-5) printf("3R %le\n", relErrorReal(5, c, zeror, zeroi));
                if (relErrorComplex(5, d, e, zeror, zeroi) > 1e-5) printf("3C %le\n", relErrorComplex(5, d, e, zeror, zeroi));
                if (relErrorComplex(5, d, -e, zeror, zeroi) > 1e-5) printf("3C %le\n", relErrorComplex(5, d, -e, zeror, zeroi));
                break;
            case 5:
                if (relErrorReal(5, a, zeror, zeroi) > 1e-5) printf("5R %le\n", relErrorReal(5, a, zeror, zeroi));
                if (relErrorReal(5, b, zeror, zeroi) > 1e-5) printf("5R %le\n", relErrorReal(5, b, zeror, zeroi));
                if (relErrorReal(5, c, zeror, zeroi) > 1e-5) printf("5R %le\n", relErrorReal(5, c, zeror, zeroi));
                if (relErrorReal(5, d, zeror, zeroi) > 1e-5) printf("5R %le\n", relErrorReal(5, d, zeror, zeroi));
                if (relErrorReal(5, e, zeror, zeroi) > 1e-5) printf("5R %le\n", relErrorReal(5, e, zeror, zeroi));
                break;
        }
    }
}

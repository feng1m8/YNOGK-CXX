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
    double coefs[3];
    double zeror[2];
    double zeroi[2];
    for (int i=0; i<100000000; i++) {
        double a = 5 - 10*drand48();
        double b = 5 - 10*drand48();
        int nReal = (lrand48() & 1) * 2;
        quadraticCoefs(0.1+drand48(), nReal, a, b, coefs);
        quadraticSolve(coefs[2], coefs[1], coefs[0], zeror, zeroi);
        switch (nReal) {
            case 0:
                if (relErrorComplex(2, a, b, zeror, zeroi) > 1e-9) printf("0C %le\n", relErrorComplex(2, a, b, zeror, zeroi));
                if (relErrorComplex(2, a, -b, zeror, zeroi) > 1e-9) printf("0C %le\n", relErrorComplex(2, a, -b, zeror, zeroi));
                break;
            case 2:
                if (relErrorReal(2, a, zeror, zeroi) > 1e-9) printf("2R %le\n", relErrorReal(2, a, zeror, zeroi));
                if (relErrorReal(2, b, zeror, zeroi) > 1e-9) printf("2R %le\n", relErrorReal(2, b, zeror, zeroi));
                break;
        }
    }
}

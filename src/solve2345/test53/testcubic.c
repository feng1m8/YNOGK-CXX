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
    double coefs[4];
    double zeror[3];
    double zeroi[3];
    for (int i=0; i<100000000; i++) {
        double a = 5 - 10*drand48();
        double b = 5 - 10*drand48();
        double c = 5 - 10*drand48();
        int nReal = 1 + (lrand48() & 1) * 2;
        cubicCoefs(0.1+drand48(), nReal, a, b, c, coefs);
        cubicSolve(coefs[3], coefs[2], coefs[1], coefs[0], zeror, zeroi);
        switch (nReal) {
            case 1:
                if (relErrorReal(3, a, zeror, zeroi) > 1e-8) printf("1R %le\n", relErrorReal(3, a, zeror, zeroi));
                if (relErrorComplex(3, b, c, zeror, zeroi) > 1e-8) printf("1C %le\n", relErrorComplex(3, b, c, zeror, zeroi));
                if (relErrorComplex(3, b, -c, zeror, zeroi) > 1e-8) printf("1C %le\n", relErrorComplex(3, b, -c, zeror, zeroi));
                break;
            case 3:
                if (relErrorReal(3, a, zeror, zeroi) > 1e-8) printf("3R %le\n", relErrorReal(3, a, zeror, zeroi));
                if (relErrorReal(3, b, zeror, zeroi) > 1e-8) printf("3R %le\n", relErrorReal(3, b, zeror, zeroi));
                if (relErrorReal(3, c, zeror, zeroi) > 1e-8) printf("3R %le\n", relErrorReal(3, c, zeror, zeroi));
                break;
        }
    }
}


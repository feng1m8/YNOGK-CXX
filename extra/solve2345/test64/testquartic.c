#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>

#include "solve2345.h"
#include "test64.h"


/**
 *
 * @author woelen
 */
    
int main(int argc, char **argv) {
    srand48(time(NULL));
    long double coefs[5];
    long double zeror[4];
    long double zeroi[4];
    for (int i=0; i<100000000; i++) {
        long double a = 5 - 10*drand48();
        long double b = 5 - 10*drand48();
        long double c = 5 - 10*drand48();
        long double d = 5 - 10*drand48();
        int nReal = (lrand48() % 3) * 2;
        quarticCoefsl(0.1+drand48(), nReal, a, b, c, d, coefs);
        quarticSolvel(coefs[4], coefs[3], coefs[2], coefs[1], coefs[0], zeror, zeroi);
        switch (nReal) {
            case 0:
                if (relErrorComplex(4, a, b, zeror, zeroi) > 1e-10) printf("0C %Le\n", relErrorComplex(4, a, b, zeror, zeroi));
                if (relErrorComplex(4, a, -b, zeror, zeroi) > 1e-10) printf("0C %Le\n", relErrorComplex(4, a, -b, zeror, zeroi));
                if (relErrorComplex(4, c, d, zeror, zeroi) > 1e-10) printf("0C %Le\n", relErrorComplex(4, c, d, zeror, zeroi));
                if (relErrorComplex(4, c, -d, zeror, zeroi) > 1e-10) printf("0C %Le\n", relErrorComplex(4, c, -d, zeror, zeroi));
                break;
            case 2:
                if (relErrorReal(4, a, zeror, zeroi) > 1e-10) printf("2R %Le\n", relErrorReal(4, a, zeror, zeroi));
                if (relErrorReal(4, b, zeror, zeroi) > 1e-10) printf("2R %Le\n", relErrorReal(4, b, zeror, zeroi));
                if (relErrorComplex(4, c, d, zeror, zeroi) > 1e-10) printf("2C %Le\n", relErrorComplex(4, c, d, zeror, zeroi));
                if (relErrorComplex(4, c, -d, zeror, zeroi) > 1e-10) printf("2C %Le\n", relErrorComplex(4, c, -d, zeror, zeroi));
                break;
            case 4:
                if (relErrorReal(4, a, zeror, zeroi) > 1e-10) printf("4R %Le\n", relErrorReal(4, a, zeror, zeroi));
                if (relErrorReal(4, b, zeror, zeroi) > 1e-10) printf("4R %Le\n", relErrorReal(4, b, zeror, zeroi));
                if (relErrorReal(4, c, zeror, zeroi) > 1e-10) printf("4R %Le\n", relErrorReal(4, c, zeror, zeroi));
                if (relErrorReal(4, d, zeror, zeroi) > 1e-10) printf("4R %Le\n", relErrorReal(4, d, zeror, zeroi));
                break;
        }
    }
}

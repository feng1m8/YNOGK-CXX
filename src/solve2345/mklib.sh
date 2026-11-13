#!/usr/bin/sh
cd src
rm -f *.o *.a

gcc -I.. -Wall -ffast-math -O3 -DFP_PREC=53 -o d2.o -c quadratic.c
gcc -I.. -Wall -ffast-math -O3 -DFP_PREC=53 -o d3.o -c cubic.c
gcc -I.. -Wall -ffast-math -O3 -DFP_PREC=53 -o d4.o -c quartic.c
gcc -I.. -Wall -ffast-math -O3 -DFP_PREC=53 -o d5.o -c quintic.c

gcc -I.. -Wall -ffast-math -O3 -DFP_PREC=64 -o ld2.o -c quadratic.c
gcc -I.. -Wall -ffast-math -O3 -DFP_PREC=64 -o ld3.o -c cubic.c
gcc -I.. -Wall -ffast-math -O3 -DFP_PREC=64 -o ld4.o -c quartic.c
gcc -I.. -Wall -ffast-math -O3 -DFP_PREC=64 -o ld5.o -c quintic.c
ar r libsolve2345.a d2.o d3.o d4.o d5.o ld2.o ld3.o ld4.o ld5.o
rm -f *.o

ranlib libsolve2345.a
mv libsolve2345.a ..

cd ..

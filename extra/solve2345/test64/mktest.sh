#!/usr/bin/sh
gcc -O3 -I.. -o testquadratic testquadratic.c ../libsolve2345.a -lm
gcc -O3 -I.. -o testcubic testcubic.c ../libsolve2345.a -lm
gcc -O3 -I.. -o testquartic testquartic.c ../libsolve2345.a -lm
gcc -O3 -I.. -o testquintic testquintic.c ../libsolve2345.a -lm

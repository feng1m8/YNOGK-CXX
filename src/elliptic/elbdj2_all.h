#ifndef ELBDJ2_ALL_H
#define ELBDJ2_ALL_H

void elbdj2(double phi, double phic, double n, double mc, double *b, double *d, double *j);
void celbdj(double nc, double mc, double *celb, double *celd, double *celj);
void elsbdj(double s, double n, double mc, double *b, double *d, double *j);
double serj(double y, double n, double m);
double uatan(double t, double h);

#endif

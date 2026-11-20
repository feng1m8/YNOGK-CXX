#ifndef YNOGK_CXX_ELLIPTIC_ELBDJ2_H
#define YNOGK_CXX_ELLIPTIC_ELBDJ2_H

#ifdef __cplusplus
extern "C"
{
#endif

    void elbdj2(double phi, double phic, double n, double mc, double *b, double *d, double *j);
    void celbdj(double nc, double mc0, double *celb, double *celd, double *celj);
    void elsbdj(double s, double n, double mc, double *b, double *d, double *j);
    double serj(double y, double n, double m);
    double uatan(double t, double h);

#ifdef __cplusplus
}
#endif

#endif

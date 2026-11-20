#ifndef YNOGK_CXX_ELLIPTIC_ELBD_H
#define YNOGK_CXX_ELLIPTIC_ELBD_H

#ifdef __cplusplus
extern "C"
{
#endif

    void elbd(double phi, double phic, double mc, double *b, double *d);
    void elsbd(double s, double mc, double *b, double *d);
    void elcbd(double c, double mc, double *b, double *d);
    void serbd(double y, double m, double *b, double *d);

#ifdef __cplusplus
}
#endif

#endif

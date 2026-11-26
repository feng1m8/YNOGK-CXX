#ifndef GELBD_H
#define GELBD_H

void gelbd(double phi, double mc, double *elb, double *eld);
void elbd(double phi, double phic, double mc, double *b, double *d);
void elsbd(double s, double mc, double *b, double *d);
void elcbd(double c, double mc, double *b, double *d);
void serbd(double y, double m, double *b, double *d);

void rgelbd(float phi, float mc, float *elb, float *eld);
void relbd(float phi, float phic, float mc, float *b, float *d);
void relsbd(float s, float mc, float *b, float *d);
void relcbd(float c, float mc, float *b, float *d);
void rserbd(float y, float m, float *b, float *d);

#endif

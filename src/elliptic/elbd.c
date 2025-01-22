#include <fenv.h>
#include <math.h>

#include "cei2.h"
#include "elbd.h"

void elbd(double phi, double phic, double mc, double *b, double *d)
{
    if (phi < 1.25)
        return elsbd(sin(phi), mc, b, d);

    double c = sin(phic);
    double x = c * c;
    double d2 = mc + (1.0 - mc) * x;

    if (x < 0.9 * d2)
    {
        double z = c / sqrt(d2);
        elsbd(z, mc, b, d);

        double sz = z * sqrt(1.0 - x);
        *b = ceib(mc) - (*b - sz);
        *d = ceid(mc) - (*d + sz);
    }
    else
    {
        double v = mc * (1.0 - x);
        if (v < x * d2)
            return elcbd(c, mc, b, d);
        else
        {
            double t2 = (1.0 - x) / d2;
            elcbd(sqrt(mc * t2), mc, b, d);

            double sz = c * sqrt(t2);
            *b = ceib(mc) - (*b - sz);
            *d = ceid(mc) - (*d + sz);
        }
    }
}

void elsbd(double s, double mc, double *b, double *d)
{
    double m = 1.0 - mc;
    double del = 0.04094 - 0.00652 * m;
    double y = s * s;
    if (y < del)
    {
        serbd(y, m, b, d);
        *b = s * *b;
        *d = s * y * *d;
        return;
    }
    double ss[11];
    double yy[11];
    ss[0] = s;
    int j;
    for (j = 1; j <= 10; ++j)
    {
        y = y / ((1.0 + sqrt(1.0 - y)) * (1.0 + sqrt(1.0 - m * y)));
        yy[j] = y;
        ss[j] = sqrt(y);
        if (y < del)
            goto label1;
    }
    feraiseexcept(FE_INVALID);
label1:

    serbd(y, m, b, d);
    *b = ss[j] * *b;
    *d = ss[j] * y * *d;
    for (int i = 1; i <= j; ++i)
    {
        double sy = ss[j - i] * yy[j - i + 1];
        *b = *b + (*b - sy);
        *d = *d + (*d + sy);
    }
}

void elcbd(double c, double mc, double *b, double *dx)
{
    double x = c * c;
    double y = 1.0 - x;
    double s = sqrt(y);

    if (x > 0.1)
        return elsbd(s, mc, b, dx);

    double m = 1.0 - mc;
    double ss[11];
    double yy[11];
    ss[0] = s;
    int j;
    for (j = 1; j <= 10; ++j)
    {
        double d = sqrt(mc + m * x);
        x = (c + d) / (1.0 + d);
        y = 1.0 - x;
        yy[j] = y;
        ss[j] = sqrt(y);
        if (x > 0.1)
            goto label1;

        c = sqrt(x);
    }
    feraiseexcept(FE_INVALID);

label1:
    elsbd(ss[j], mc, b, dx);

    for (int i = 1; i <= j; ++i)
    {
        double sy = ss[j - i] * yy[j - i + 1];
        *b = *b + (*b - sy);
        *dx = *dx + (*dx + sy);
    }
}

void serbd(double y, double m, double *b, double *d)
{
    const double F10 = 1.0 / 6.0;
    const double F20 = 3.0 / 40.0;
    const double F21 = 2.0 / 40.0;
    const double F30 = 5.0 / 112.0;
    const double F31 = 3.0 / 112.0;
    const double F40 = 35.0 / 1152.0;
    const double F41 = 20.0 / 1152.0;
    const double F42 = 18.0 / 1152.0;
    const double F50 = 63.0 / 2816.0;
    const double F51 = 35.0 / 2816.0;
    const double F52 = 30.0 / 2816.0;
    const double F60 = 231.0 / 13312.0;
    const double F61 = 126.0 / 13312.0;
    const double F62 = 105.0 / 13312.0;
    const double F63 = 100.0 / 13312.0;
    const double F70 = 429.0 / 30720.0;
    const double F71 = 231.0 / 30720.0;
    const double F72 = 189.0 / 30720.0;
    const double F73 = 175.0 / 30720.0;
    const double F80 = 6435.0 / 557056.0;
    const double F81 = 3432.0 / 557056.0;
    const double F82 = 2722.0 / 557056.0;
    const double F83 = 2520.0 / 557056.0;
    const double F84 = 2450.0 / 557056.0;
    const double F90 = 12155.0 / 1245184.0;
    const double F91 = 6435.0 / 1245184.0;
    const double F92 = 5148.0 / 1245184.0;
    const double F93 = 4620.0 / 1245184.0;
    const double F94 = 4410.0 / 1245184.0;
    const double FA0 = 46189.0 / 5505024.0;
    const double FA1 = 24310.0 / 5505024.0;
    const double FA2 = 19305.0 / 5505024.0;
    const double FA3 = 17160.0 / 5505024.0;
    const double FA4 = 16170.0 / 5505024.0;
    const double FA5 = 15876.0 / 5505024.0;
    const double FB0 = 88179.0 / 12058624.0;
    const double FB1 = 46189.0 / 12058624.0;
    const double FB2 = 36465.0 / 12058624.0;
    const double FB3 = 32175.0 / 12058624.0;
    const double FB4 = 30030.0 / 12058624.0;
    const double FB5 = 29106.0 / 12058624.0;

    const double A1 = 3.0 / 5.0;
    const double A2 = 5.0 / 7.0;
    const double A3 = 7.0 / 9.0;
    const double A4 = 9.0 / 11.0;
    const double A5 = 11.0 / 13.0;
    const double A6 = 13.0 / 15.0;
    const double A7 = 15.0 / 17.0;
    const double A8 = 17.0 / 19.0;
    const double A9 = 19.0 / 21.0;
    const double AA = 21.0 / 23.0;
    const double AB = 23.0 / 25.0;

    double F1 = F10 + m * F10;
    double F2 = F20 + m * (F21 + m * F20);
    double F3 = F30 + m * (F31 + m * (F31 + m * F30));
    double F4 = F40 + m * (F41 + m * (F42 + m * (F41 + m * F40)));
    double F5 = F50 + m * (F51 + m * (F52 + m * (F52 + m * (F51 + m * F50))));
    double F6 = F60 + m * (F61 + m * (F62 + m * (F63 + m * (F62 + m * (F61 + m * F60)))));
    double F7 = F70 + m * (F71 + m * (F72 + m * (F73 + m * (F73 + m * (F72 + m * (F71 + m * F70))))));
    double F8 = F80 + m * (F81 + m * (F82 + m * (F83 + m * (F84 + m * (F83 + m * (F82 + m * (F81 + m * F80)))))));
    double F9 = F90 + m * (F91 + m * (F92 + m * (F93 + m * (F94 + m * (F94 + m * (F93 + m * (F92 + m * (F91 + m * F90))))))));
    double FA = FA0 + m * (FA1 + m * (FA2 + m * (FA3 + m * (FA4 + m * (FA5 + m * (FA4 + m * (FA3 + m * (FA2 + m * (FA1 + m * FA0)))))))));
    double FB = FB0 + m * (FB1 + m * (FB2 + m * (FB3 + m * (FB4 + m * (FB5 + m * (FB5 + m * (FB4 + m * (FB3 + m * (FB2 + m * (FB1 + m * FB0))))))))));

    const double D0 = 1.0 / 3.0;
    double D1 = F1 * A1;
    double D2 = F2 * A2;
    double D3 = F3 * A3;
    double D4 = F4 * A4;
    double D5 = F5 * A5;
    double D6 = F6 * A6;
    double D7 = F7 * A7;
    double D8 = F8 * A8;
    double D9 = F9 * A9;
    double DA = FA * AA;
    double DB = FB * AB;
    *d = D0 + y * (D1 + y * (D2 + y * (D3 + y * (D4 + y * (D5 + y * (D6 + y * (D7 + y * (D8 + y * (D9 + y * (DA + y * DB))))))))));

    double B1 = F1 - D0;
    double B2 = F2 - D1;
    double B3 = F3 - D2;
    double B4 = F4 - D3;
    double B5 = F5 - D4;
    double B6 = F6 - D5;
    double B7 = F7 - D6;
    double B8 = F8 - D7;
    double B9 = F9 - D8;
    double BA = FA - D9;
    double BB = FB - DA;
    *b = 1.0 + y * (B1 + y * (B2 + y * (B3 + y * (B4 + y * (B5 + y * (B6 + y * (B7 + y * (B8 + y * (B9 + y * (BA + y * BB))))))))));
}

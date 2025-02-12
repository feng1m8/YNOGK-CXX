#include <fenv.h>
#include <math.h>

#include "cei2.h"
#include "elf.h"

double elf(double phi, double phic, double mc)
{
    double m = 1.0 - mc;
    if (phi < 1.25)
        return asn(sin(phi), m);

    double c = sin(phic);
    double yc = c * c;
    double d2 = mc + m * yc;

    if (yc < 0.9 * d2)
        return ceik(mc) - asn(c / sqrt(d2), m);

    double v = mc * (1.0 - yc);
    if (v < yc * d2)
        return acn(c, mc);

    return ceik(mc) - asn(sqrt(1.0 - v / d2), m);
}

double acn(double c, double mc)
{
    double m = 1.0 - mc;
    double yc = c * c;

    if (yc > 0.5)
        return asn(sqrt(1.0 - yc), m);

    double s;
    double f = 1.0;
    for (int i = 1; i <= 20; ++i)
    {
        double d = sqrt(mc + m * yc);
        yc = (c + d) / (1.0 + d);
        f = f * 2.0;
        if (yc > 0.5)
        {
            s = sqrt(1.0 - yc);
            goto label1;
        }
        c = sqrt(yc);
    }
    feraiseexcept(FE_INVALID);

label1:
    return f * asn(s, m);
}

double asn(double s, double m)
{
    double del = 0.04094 - 0.00652 * m;
    double y = s * s;
    if (y < del)
        return s * serf(y, m);

    double f = 1.0;
    for (int j = 1; j <= 20; ++j)
    {
        y = y / ((1.0 + sqrt(1.0 - y)) * (1.0 + sqrt(1.0 - m * y)));
        f = f * 2.0;
        if (y < del)
        {
            s = sqrt(y);
            goto label1;
        }
    }
    feraiseexcept(FE_INVALID);

label1:
    return f * s * serf(y, m);
}

double serf(double y, double m)
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

    double F1 = F10 + m * F10;
    double F2 = F20 + m * (F21 + m * F20);
    double F3 = F30 + m * (F31 + m * (F31 + m * F30));
    double F4 = F40 + m * (F41 + m * (F42 + m * (F41 + m * F40)));
    double F5 = F50 + m * (F51 + m * (F52 + m * (F52 + m * (F51 + m * F50))));
    double F6 = F60 + m * (F61 + m * (F62 + m * (F63 + m * (F62 + m * (F61 + m * F60)))));
    double F7 = F70 + m * (F71 + m * (F72 + m * (F73 + m * (F73 + m * (F72 + m * (F71 + m * F70))))));
    double F8 = F80 + m * (F81 + m * (F82 + m * (F83 + m * (F84 + m * (F83 + m * (F82 + m * (F81 + m * F80)))))));
    double F9 = F90 + m * (F91 + m * (F92 + m * (F93 + m * (F94 + m * (F94 + m * (F93 + m * (F92 + m * (F91 + m * F90))))))));
    return 1.0 + y * (F1 + y * (F2 + y * (F3 + y * (F4 + y * (F5 + y * (F6 + y * (F7 + y * (F8 + y * F9))))))));
}

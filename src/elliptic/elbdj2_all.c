#include <fenv.h>
#include <iso646.h>
#include <math.h>

#include "cei2.h"
#include "elbdj2_all.h"
#include "gelbd.h"

void elbdj2(double phi, double phic, double n, double mc, double *b, double *d, double *j)
{
    if (phi < 1.249)
        return elsbdj(sin(phi), n, mc, b, d, j);

    double m = 1.0 - mc;
    double nc = 1.0 - n;
    double h = n * nc * (n - m);
    double c = sin(phic);
    double x = c * c;
    double d2 = mc + m * x;

    if (x < 0.9 * d2)
    {
        double z = c / sqrt(d2);
        elsbdj(z, n, mc, b, d, j);

        double bc, dc, jc;
        celbdj(nc, mc, &bc, &dc, &jc);

        double sz = z * sqrt(1.0 - x);
        double t = sz / nc;
        *b = bc - (*b - sz);
        *d = dc - (*d + sz);
        *j = jc - (*j + uatan(t, h));
    }
    else
    {
        double v = mc * (1.0 - x);
        if (v < x * d2)
            return elsbdj(sqrt(1.0 - x), n, mc, b, d, j);

        double t2 = (1.0 - x) / d2;
        elsbdj(sqrt(1.0 - mc * t2), n, mc, b, d, j);

        double bc, dc, jc;
        celbdj(nc, mc, &bc, &dc, &jc);

        double sz = c * sqrt(t2);
        double t = sz / nc;
        *b = bc - (*b - sz);
        *d = dc - (*d + sz);
        *j = jc - (*j + uatan(t, h));
    }
}

void celbdj(double nc, double mc, double *celb, double *celd, double *celj)
{
    static const double PI = 3.1415926535897932;

    enum
    {
        IMAX = 40,
    };

    static const double one_minus_B[] = {
        1.0 - 1.0 / (2.0 * 1.0 + 1.0),
        1.0 - 1.0 / (2.0 * 2.0 + 1.0),
        1.0 - 1.0 / (2.0 * 3.0 + 1.0),
        1.0 - 1.0 / (2.0 * 4.0 + 1.0),
        1.0 - 1.0 / (2.0 * 5.0 + 1.0),
        1.0 - 1.0 / (2.0 * 6.0 + 1.0),
        1.0 - 1.0 / (2.0 * 7.0 + 1.0),
        1.0 - 1.0 / (2.0 * 8.0 + 1.0),
        1.0 - 1.0 / (2.0 * 9.0 + 1.0),
        1.0 - 1.0 / (2.0 * 10.0 + 1.0),
        1.0 - 1.0 / (2.0 * 11.0 + 1.0),
        1.0 - 1.0 / (2.0 * 12.0 + 1.0),
        1.0 - 1.0 / (2.0 * 13.0 + 1.0),
        1.0 - 1.0 / (2.0 * 14.0 + 1.0),
        1.0 - 1.0 / (2.0 * 15.0 + 1.0),
        1.0 - 1.0 / (2.0 * 16.0 + 1.0),
        1.0 - 1.0 / (2.0 * 17.0 + 1.0),
        1.0 - 1.0 / (2.0 * 18.0 + 1.0),
        1.0 - 1.0 / (2.0 * 19.0 + 1.0),
        1.0 - 1.0 / (2.0 * 20.0 + 1.0),
        1.0 - 1.0 / (2.0 * 21.0 + 1.0),
        1.0 - 1.0 / (2.0 * 22.0 + 1.0),
        1.0 - 1.0 / (2.0 * 23.0 + 1.0),
        1.0 - 1.0 / (2.0 * 24.0 + 1.0),
        1.0 - 1.0 / (2.0 * 25.0 + 1.0),
        1.0 - 1.0 / (2.0 * 26.0 + 1.0),
        1.0 - 1.0 / (2.0 * 27.0 + 1.0),
        1.0 - 1.0 / (2.0 * 28.0 + 1.0),
        1.0 - 1.0 / (2.0 * 29.0 + 1.0),
        1.0 - 1.0 / (2.0 * 30.0 + 1.0),
        1.0 - 1.0 / (2.0 * 31.0 + 1.0),
        1.0 - 1.0 / (2.0 * 32.0 + 1.0),
        1.0 - 1.0 / (2.0 * 33.0 + 1.0),
        1.0 - 1.0 / (2.0 * 34.0 + 1.0),
        1.0 - 1.0 / (2.0 * 35.0 + 1.0),
        1.0 - 1.0 / (2.0 * 36.0 + 1.0),
        1.0 - 1.0 / (2.0 * 37.0 + 1.0),
        1.0 - 1.0 / (2.0 * 38.0 + 1.0),
        1.0 - 1.0 / (2.0 * 39.0 + 1.0),
        1.0 - 1.0 / (2.0 * 40.0 + 1.0),
    };

    static const double two_B_minus_1[] = {
        -1.0 + 2.0 / (2.0 * 1.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 2.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 3.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 4.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 5.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 6.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 7.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 8.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 9.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 10.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 11.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 12.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 13.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 14.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 15.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 16.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 17.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 18.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 19.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 20.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 21.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 22.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 23.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 24.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 25.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 26.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 27.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 28.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 29.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 30.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 31.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 32.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 33.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 34.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 35.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 36.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 37.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 38.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 39.0 + 1.0),
        -1.0 + 2.0 / (2.0 * 40.0 + 1.0),
    };

    if (isnan(nc))
    {
        *celb = ceib(mc);
        *celd = ceid(mc);
        *celj = nc;
        return;
    }
    else if (mc > 1.0)
    {
        if (isinf(mc))
        {
            *celb = *celd = *celj = 0.0;
            return;
        }

        mc = 1.0 / mc;
        double kc = sqrt(mc);
        nc = nc * mc;

        celbdj(nc, mc, celd, celb, celj);
        *celb = *celb * kc;
        *celd = *celd * kc;
        *celj = *celj * kc * mc;
        return;
    }
    else if (mc == 1.0)
    {
        *celb = *celd = 0.25 * PI;

        if (nc >= 0.0)
            *celj = 0.5 * PI / (nc + sqrt(nc));
        else
            *celj = 0.5 * PI / (nc - 1.0);

        return;
    }
    else if (mc > 0.0)
    {
        if (isinf(nc))
        {
            *celb = ceib(mc);
            *celd = ceid(mc);
            *celj = 0.0;
            return;
        }
        else if (nc > mc)
        {
            double m = 1.0 - mc;
            double m1 = 1.0 + m;
            double n = 1.0 - nc;

            *celb = ceib(mc);
            *celd = ceid(mc);
            double celk = *celb + *celd;

            double c[IMAX];
            double d[IMAX];
            double y[IMAX + 1];
            y[0] = n / m;

            int is = 0;
            if (y[0] > 0.5)
            {
                double x = 1.0 - y[0];
                for (int i = 0; i <= IMAX - 1; ++i)
                {
                    c[i] = sqrt(x);
                    d[i] = sqrt(mc + m * x);
                    x = (c[i] + d[i]) / (1.0 + d[i]);

                    if (x > 0.5)
                    {
                        y[i + 1] = y[i] / ((1.0 + c[i]) * (1.0 + d[i]));
                        is = i + 1;
                        break;
                    }

                    y[i + 1] = 1.0 - x;
                }

                if (is == 0)
                {
                    feraiseexcept(FE_INVALID);
                    is = IMAX;
                }
            }

            for (int i = is; i <= IMAX - 1; ++i)
            {
                c[i] = sqrt(1.0 - y[i]);
                d[i] = sqrt(1.0 - m * y[i]);
                y[i + 1] = y[i] / ((1.0 + c[i]) * (1.0 + d[i]));

                if (fabs(y[i + 1]) < 0.325)
                {
                    is = i;
                    break;
                }
            }

            if (is >= IMAX)
            {
                feraiseexcept(FE_INVALID);
                is = IMAX - 1;
            }

            double a[IMAX + 2];
            a[0] = celk / m;
            a[1] = *celd;

            *celj = a[1];
            double yi = 1.0;
            double ye = y[is + 1];

            for (int i = 1; i <= IMAX; ++i)
            {
                yi = yi * ye;
                a[i + 1] = one_minus_B[i - 1] * m1 * a[i] + two_B_minus_1[i - 1] * m * a[i - 1];
                double dj = a[i + 1] * yi;

                if (fabs(dj) < 1.11e-16 * fabs(*celj))
                    break;

                *celj = *celj + dj;
            }

            for (int i = is; i >= 0; --i)
            {
                *celj = (2.0 * (c[i] + d[i]) * *celj - y[i] * celk) / (c[i] * d[i] * (1.0 + c[i]) * (1.0 + d[i]));
            }

            return;
        }
        else if (nc > 0.0)
        {
            nc = 1.0 / nc;
            double n2c = mc * nc;

            celbdj(n2c, mc, celb, celd, celj);
            double celk = *celb + *celd;
            *celj = (celk - n2c * *celj) * nc;
            return;
        }
        else if (nc < 0.0)
        {
            double n = 1.0 - nc;
            n = 1.0 / n;
            double m = 1.0 - mc;
            double n1 = m * n;
            nc = 1.0 - n1;

            celbdj(nc, mc, celb, celd, celj);
            double celk = *celb + *celd;
            *celj = -(celk + n1 * *celj) * n;
            return;
        }
        else
        {
            *celb = ceib(mc);
            *celd = ceid(mc);
            *celj = HUGE_VAL;
            return;
        }
    }
    else if (mc == 0.0)
    {
        *celb = 1.0;
        *celd = *celj = HUGE_VAL;
        return;
    }
    else
    {
        if (not isnan(mc))
            feraiseexcept(FE_INVALID);

        *celb = *celd = *celj = nan("");
        return;
    }
}

void elsbdj(double s, double n, double mc, double *b, double *d, double *j)
{
    static const double del = 0.01622;

    double m = 1.0 - mc;
    double y = s * s;

    if (y < del)
    {
        serbd(y, m, b, d);
        *b = s * *b;
        *d = s * *d * y;
        *j = s * serj(y, n, m);
        return;
    }

    double yy[11];
    double ss[11];
    double cd[10];
    yy[0] = y;
    ss[0] = s;

    int i;
    for (i = 1; i <= 10; ++i)
    {
        double c = sqrt(1.0 - y);
        double d = sqrt(1.0 - m * y);
        y = y / ((1.0 + c) * (1.0 + d));

        yy[i] = y;
        ss[i] = sqrt(y);
        cd[i - 1] = c * d;

        if (y < del)
            break;
    }

    if (i > 10)
    {
        feraiseexcept(FE_INVALID);
        i = 10;
    }

    serbd(y, m, b, d);
    *b = ss[i] * *b;
    *d = ss[i] * *d * y;
    *j = ss[i] * serj(y, n, m);

    double h = n * (1.0 - n) * (n - m);

    for (int k = i - 1; k >= 0; --k)
    {
        double sy = ss[k] * yy[k + 1];
        double t = sy / (1.0 - n * (yy[k] - yy[k + 1] * cd[k]));

        *b = 2.0 * *b - sy;
        *d = 2.0 * *d + sy;
        *j = 2.0 * *j + uatan(t, h);
    }
}

double serj(double y, double n, double m)
{
    const double J100 = 1.0 / 3.0;

    const double J200 = 1.0 / 10.0;
    const double J201 = 2.0 / 10.0;
    const double J210 = 1.0 / 10.0;

    const double J300 = 3.0 / 56.0;
    const double J301 = 4.0 / 56.0;
    const double J302 = 8.0 / 56.0;
    const double J310 = 2.0 / 56.0;
    const double J311 = 4.0 / 56.0;
    const double J320 = 3.0 / 56.0;

    const double J400 = 5.0 / 144.0;
    const double J401 = 6.0 / 144.0;
    const double J402 = 8.0 / 144.0;
    const double J403 = 16.0 / 144.0;
    const double J410 = 3.0 / 144.0;
    const double J411 = 4.0 / 144.0;
    const double J412 = 8.0 / 144.0;
    const double J420 = 3.0 / 144.0;
    const double J421 = 6.0 / 144.0;
    const double J430 = 5.0 / 144.0;

    const double J500 = 35.0 / 1408.0;
    const double J501 = 40.0 / 1408.0;
    const double J502 = 48.0 / 1408.0;
    const double J503 = 64.0 / 1408.0;
    const double J504 = 128.0 / 1408.0;
    const double J510 = 20.0 / 1408.0;
    const double J511 = 24.0 / 1408.0;
    const double J512 = 32.0 / 1408.0;
    const double J513 = 64.0 / 1408.0;
    const double J520 = 18.0 / 1408.0;
    const double J521 = 24.0 / 1408.0;
    const double J522 = 48.0 / 1408.0;
    const double J530 = 20.0 / 1408.0;
    const double J531 = 40.0 / 1408.0;
    const double J540 = 35.0 / 1408.0;

    const double J600 = 63.0 / 3328.0;
    const double J601 = 70.0 / 3328.0;
    const double J602 = 80.0 / 3328.0;
    const double J603 = 96.0 / 3328.0;
    const double J604 = 128.0 / 3328.0;
    const double J605 = 256.0 / 3328.0;
    const double J610 = 35.0 / 3328.0;
    const double J611 = 40.0 / 3328.0;
    const double J612 = 48.0 / 3328.0;
    const double J613 = 64.0 / 3328.0;
    const double J614 = 128.0 / 3328.0;
    const double J620 = 30.0 / 3328.0;
    const double J621 = 36.0 / 3328.0;
    const double J622 = 48.0 / 3328.0;
    const double J623 = 96.0 / 3328.0;
    const double J630 = 30.0 / 3328.0;
    const double J631 = 40.0 / 3328.0;
    const double J632 = 80.0 / 3328.0;
    const double J640 = 35.0 / 3328.0;
    const double J641 = 70.0 / 3328.0;
    const double J650 = 63.0 / 3328.0;

    const double J700 = 231.0 / 15360.0;
    const double J701 = 252.0 / 15360.0;
    const double J702 = 280.0 / 15360.0;
    const double J703 = 320.0 / 15360.0;
    const double J704 = 384.0 / 15360.0;
    const double J705 = 512.0 / 15360.0;
    const double J706 = 1024.0 / 15360.0;
    const double J710 = 126.0 / 15360.0;
    const double J711 = 140.0 / 15360.0;
    const double J712 = 160.0 / 15360.0;
    const double J713 = 192.0 / 15360.0;
    const double J714 = 256.0 / 15360.0;
    const double J715 = 512.0 / 15360.0;
    const double J720 = 105.0 / 15360.0;
    const double J721 = 120.0 / 15360.0;
    const double J722 = 144.0 / 15360.0;
    const double J723 = 192.0 / 15360.0;
    const double J724 = 384.0 / 15360.0;
    const double J730 = 100.0 / 15360.0;
    const double J731 = 120.0 / 15360.0;
    const double J732 = 160.0 / 15360.0;
    const double J733 = 320.0 / 15360.0;
    const double J740 = 105.0 / 15360.0;
    const double J741 = 140.0 / 15360.0;
    const double J742 = 280.0 / 15360.0;
    const double J750 = 126.0 / 15360.0;
    const double J751 = 252.0 / 15360.0;
    const double J760 = 231.0 / 15360.0;

    const double J800 = 429.0 / 34816.0;
    const double J801 = 462.0 / 34816.0;
    const double J802 = 504.0 / 34816.0;
    const double J803 = 560.0 / 34816.0;
    const double J804 = 640.0 / 34816.0;
    const double J805 = 768.0 / 34816.0;
    const double J806 = 1024.0 / 34816.0;
    const double J807 = 2048.0 / 34816.0;
    const double J810 = 231.0 / 34816.0;
    const double J811 = 252.0 / 34816.0;
    const double J812 = 280.0 / 34816.0;
    const double J813 = 320.0 / 34816.0;
    const double J814 = 384.0 / 34816.0;
    const double J815 = 512.0 / 34816.0;
    const double J816 = 1024.0 / 34816.0;
    const double J820 = 189.0 / 34816.0;
    const double J821 = 210.0 / 34816.0;
    const double J822 = 240.0 / 34816.0;
    const double J823 = 288.0 / 34816.0;
    const double J824 = 284.0 / 34816.0;
    const double J825 = 768.0 / 34816.0;
    const double J830 = 175.0 / 34816.0;
    const double J831 = 200.0 / 34816.0;
    const double J832 = 240.0 / 34816.0;
    const double J833 = 320.0 / 34816.0;
    const double J834 = 640.0 / 34816.0;
    const double J840 = 175.0 / 34816.0;
    const double J841 = 210.0 / 34816.0;
    const double J842 = 280.0 / 34816.0;
    const double J843 = 560.0 / 34816.0;
    const double J850 = 189.0 / 34816.0;
    const double J851 = 252.0 / 34816.0;
    const double J852 = 504.0 / 34816.0;
    const double J860 = 231.0 / 34816.0;
    const double J861 = 462.0 / 34816.0;
    const double J870 = 429.0 / 34816.0;

    const double J900 = 6435.0 / 622592.0;
    const double J901 = 6864.0 / 622592.0;
    const double J902 = 7392.0 / 622592.0;
    const double J903 = 8064.0 / 622592.0;
    const double J904 = 8960.0 / 622592.0;
    const double J905 = 10240.0 / 622592.0;
    const double J906 = 12288.0 / 622592.0;
    const double J907 = 16384.0 / 622592.0;
    const double J908 = 32768.0 / 622592.0;
    const double J910 = 3432.0 / 622592.0;
    const double J911 = 3696.0 / 622592.0;
    const double J912 = 4032.0 / 622592.0;
    const double J913 = 4480.0 / 622592.0;
    const double J914 = 5120.0 / 622592.0;
    const double J915 = 6144.0 / 622592.0;
    const double J916 = 8192.0 / 622592.0;
    const double J917 = 16384.0 / 622592.0;
    const double J920 = 2772.0 / 622592.0;
    const double J921 = 3024.0 / 622592.0;
    const double J922 = 3360.0 / 622592.0;
    const double J923 = 3840.0 / 622592.0;
    const double J924 = 4608.0 / 622592.0;
    const double J925 = 6144.0 / 622592.0;
    const double J926 = 12288.0 / 622592.0;
    const double J930 = 2520.0 / 622592.0;
    const double J931 = 2800.0 / 622592.0;
    const double J932 = 3200.0 / 622592.0;
    const double J933 = 3840.0 / 622592.0;
    const double J934 = 5120.0 / 622592.0;
    const double J935 = 10240.0 / 622592.0;
    const double J940 = 2450.0 / 622592.0;
    const double J941 = 2800.0 / 622592.0;
    const double J942 = 3360.0 / 622592.0;
    const double J943 = 4480.0 / 622592.0;
    const double J944 = 8960.0 / 622592.0;
    const double J950 = 2520.0 / 622592.0;
    const double J951 = 3024.0 / 622592.0;
    const double J952 = 4032.0 / 622592.0;
    const double J953 = 8064.0 / 622592.0;
    const double J960 = 2772.0 / 622592.0;
    const double J961 = 3696.0 / 622592.0;
    const double J962 = 7392.0 / 622592.0;
    const double J970 = 3432.0 / 622592.0;
    const double J971 = 6864.0 / 622592.0;
    const double J980 = 6435.0 / 622592.0;

    const double JA00 = 12155.0 / 1376256.0;
    const double JA01 = 12870.0 / 1376256.0;
    const double JA02 = 13728.0 / 1376256.0;
    const double JA03 = 14784.0 / 1376256.0;
    const double JA04 = 16128.0 / 1376256.0;
    const double JA05 = 17920.0 / 1376256.0;
    const double JA06 = 20480.0 / 1376256.0;
    const double JA07 = 24576.0 / 1376256.0;
    const double JA08 = 32768.0 / 1376256.0;
    const double JA09 = 65536.0 / 1376256.0;
    const double JA10 = 6435.0 / 1376256.0;
    const double JA11 = 6864.0 / 1376256.0;
    const double JA12 = 7392.0 / 1376256.0;
    const double JA13 = 8064.0 / 1376256.0;
    const double JA14 = 8960.0 / 1376256.0;
    const double JA15 = 10240.0 / 1376256.0;
    const double JA16 = 12288.0 / 1376256.0;
    const double JA17 = 16384.0 / 1376256.0;
    const double JA18 = 32768.0 / 1376256.0;
    const double JA20 = 5148.0 / 1376256.0;
    const double JA21 = 5544.0 / 1376256.0;
    const double JA22 = 6048.0 / 1376256.0;
    const double JA23 = 6720.0 / 1376256.0;
    const double JA24 = 7680.0 / 1376256.0;
    const double JA25 = 9216.0 / 1376256.0;
    const double JA26 = 12288.0 / 1376256.0;
    const double JA27 = 24576.0 / 1376256.0;
    const double JA30 = 4620.0 / 1376256.0;
    const double JA31 = 5040.0 / 1376256.0;
    const double JA32 = 5600.0 / 1376256.0;
    const double JA33 = 6400.0 / 1376256.0;
    const double JA34 = 7680.0 / 1376256.0;
    const double JA35 = 10240.0 / 1376256.0;
    const double JA36 = 20480.0 / 1376256.0;
    const double JA40 = 4410.0 / 1376256.0;
    const double JA41 = 4900.0 / 1376256.0;
    const double JA42 = 5600.0 / 1376256.0;
    const double JA43 = 6720.0 / 1376256.0;
    const double JA44 = 8960.0 / 1376256.0;
    const double JA45 = 17920.0 / 1376256.0;
    const double JA50 = 4410.0 / 1376256.0;
    const double JA51 = 5040.0 / 1376256.0;
    const double JA52 = 6048.0 / 1376256.0;
    const double JA53 = 8064.0 / 1376256.0;
    const double JA54 = 16128.0 / 1376256.0;
    const double JA60 = 4620.0 / 1376256.0;
    const double JA61 = 5544.0 / 1376256.0;
    const double JA62 = 7392.0 / 1376256.0;
    const double JA63 = 14784.0 / 1376256.0;
    const double JA70 = 5148.0 / 1376256.0;
    const double JA71 = 6864.0 / 1376256.0;
    const double JA72 = 13728.0 / 1376256.0;
    const double JA80 = 6435.0 / 1376256.0;
    const double JA81 = 12870.0 / 1376256.0;
    const double JA90 = 12155.0 / 1376256.0;

    const double J1 = J100;
    double J2 = J200 + n * J201 + m * J210;
    double J3 = J300 + n * (J301 + n * J302) + m * (J310 + n * J311 + m * J320);
    double J4 = J400 + n * (J401 + n * (J402 + n * J403)) + m * (J410 + n * (J411 + n * J412) + m * (J420 + n * J421 + m * J430));
    double J5 = J500 + n * (J501 + n * (J502 + n * (J503 + n * J504))) + m * (J510 + n * (J511 + n * (J512 + n * J513)) + m * (J520 + n * (J521 + n * J522) + m * (J530 + n * J531 + m * J540)));
    if (y <= 6.0369310e-4)
        return y * (J1 + y * (J2 + y * (J3 + y * (J4 + y * J5))));

    double J6 = J600 + n * (J601 + n * (J602 + n * (J603 + n * (J604 + n * J605)))) + m * (J610 + n * (J611 + n * (J612 + n * (J613 + n * J614))) + m * (J620 + n * (J621 + n * (J622 + n * J623)) + m * (J630 + n * (J631 + n * J632) + m * (J640 + n * J641 + m * J650))));
    if (y <= 2.0727505e-3)
        return y * (J1 + y * (J2 + y * (J3 + y * (J4 + y * (J5 + y * J6)))));

    double J7 = J700 + n * (J701 + n * (J702 + n * (J703 + n * (J704 + n * (J705 + n * J706))))) + m * (J710 + n * (J711 + n * (J712 + n * (J713 + n * (J714 + n * J715)))) + m * (J720 + n * (J721 + n * (J722 + n * (J723 + n * J724))) + m * (J730 + n * (J731 + n * (J732 + n * J733)) + m * (J740 + n * (J741 + n * J742) + m * (J750 + n * J751 + m * J760)))));
    if (y <= 5.0047026e-3)
        return y * (J1 + y * (J2 + y * (J3 + y * (J4 + y * (J5 + y * (J6 + y * J7))))));

    double J8 = J800 + n * (J801 + n * (J802 + n * (J803 + n * (J804 + n * (J805 + n * (J806 + n * J807)))))) + m * (J810 + n * (J811 + n * (J812 + n * (J813 + n * (J814 + n * (J815 + n * J816))))) + m * (J820 + n * (J821 + n * (J822 + n * (J823 + n * (J824 + n * J825)))) + m * (J830 + n * (J831 + n * (J832 + n * (J833 + n * J834))) + m * (J840 + n * (J841 + n * (J842 + n * J843)) + m * (J850 + n * (J851 + n * J852) + m * (J860 + n * J861 + m * J870))))));
    if (y <= 9.6961652e-3)
        return y * (J1 + y * (J2 + y * (J3 + y * (J4 + y * (J5 + y * (J6 + y * (J7 + y * J8)))))));

    double J9 = J900 + n * (J901 + n * (J902 + n * (J903 + n * (J904 + n * (J905 + n * (J906 + n * (J907 + n * J908))))))) + m * (J910 + n * (J911 + n * (J912 + n * (J913 + n * (J914 + n * (J915 + n * (J916 + n * J917)))))) + m * (J920 + n * (J921 + n * (J922 + n * (J923 + n * (J924 + n * (J925 + n * J926))))) + m * (J930 + n * (J931 + n * (J932 + n * (J933 + n * (J934 + n * J935)))) + m * (J940 + n * (J941 + n * (J942 + n * (J943 + n * J944))) + m * (J950 + n * (J951 + n * (J952 + n * J953)) + m * (J960 + n * (J961 + n * J962) + m * (J970 + n * J971 + m * J980)))))));
    if (y <= 1.6220210e-2)
        return y * (J1 + y * (J2 + y * (J3 + y * (J4 + y * (J5 + y * (J6 + y * (J7 + y * (J8 + y * J9))))))));

    double JA = JA00 + n * (JA01 + n * (JA02 + n * (JA03 + n * (JA04 + n * (JA05 + n * (JA06 + n * (JA07 + n * (JA08 + n * JA09)))))))) + m * (JA10 + n * (JA11 + n * (JA12 + n * (JA13 + n * (JA14 + n * (JA15 + n * (JA16 + n * (JA17 + n * JA18))))))) + m * (JA20 + n * (JA21 + n * (JA22 + n * (JA23 + n * (JA24 + n * (JA25 + n * (JA26 + n * JA27)))))) + m * (JA30 + n * (JA31 + n * (JA32 + n * (JA33 + n * (JA34 + n * (JA35 + n * JA36))))) + m * (JA40 + n * (JA41 + n * (JA42 + n * (JA43 + n * (JA44 + n * JA45)))) + m * (JA50 + n * (JA51 + n * (JA52 + n * (JA53 + n * JA54))) + m * (JA60 + n * (JA61 + n * (JA62 + n * JA63)) + m * (JA70 + n * (JA71 + n * JA72) + m * (JA80 + n * JA81 + m * JA90))))))));
    return y * (J1 + y * (J2 + y * (J3 + y * (J4 + y * (J5 + y * (J6 + y * (J7 + y * (J8 + y * (J9 + y * JA)))))))));
}

double uatan(double t, double h)
{
    static const double A3 = 1.0 / 3.0;
    static const double A5 = 1.0 / 5.0;
    static const double A7 = 1.0 / 7.0;
    static const double A9 = 1.0 / 9.0;
    static const double A11 = 1.0 / 11.0;
    static const double A13 = 1.0 / 13.0;
    static const double A15 = 1.0 / 15.0;
    static const double A17 = 1.0 / 17.0;
    static const double A19 = 1.0 / 19.0;
    static const double A21 = 1.0 / 21.0;
    static const double A23 = 1.0 / 23.0;
    static const double A25 = 1.0 / 25.0;

    double z = -h * t * t;
    double a = fabs(z);

    if (a < 3.3306691e-16)
        return t;
    else if (a < 2.3560805e-8)
        return t * (1.0 + z * A3);
    else if (a < 9.1939631e-6)
        return t * (1.0 + z * (A3 + z * A5));
    else if (a < 1.7779240e-4)
        return t * (1.0 + z * (A3 + z * (A5 + z * A7)));
    else if (a < 1.0407839e-3)
        return t * (1.0 + z * (A3 + z * (A5 + z * (A7 + z * A9))));
    else if (a < 3.3616998e-3)
        return t * (1.0 + z * (A3 + z * (A5 + z * (A7 + z * (A9 + z * A11)))));
    else if (a < 7.7408014e-3)
        return t * (1.0 + z * (A3 + z * (A5 + z * (A7 + z * (A9 + z * (A11 + z * A13))))));
    else if (a < 1.4437181e-2)
        return t * (1.0 + z * (A3 + z * (A5 + z * (A7 + z * (A9 + z * (A11 + z * (A13 + z * A15)))))));
    else if (a < 2.3407312e-2)
        return t * (1.0 + z * (A3 + z * (A5 + z * (A7 + z * (A9 + z * (A11 + z * (A13 + z * (A15 + z * A17))))))));
    else if (a < 3.4416203e-2)
        return t * (1.0 + z * (A3 + z * (A5 + z * (A7 + z * (A9 + z * (A11 + z * (A13 + z * (A15 + z * (A17 + z * A19)))))))));
    else if (z < 0.0)
    {
        double r = sqrt(h);
        return atan(r * t) / r;
    }
    else if (a < 4.7138547e-2)
        return t * (1.0 + z * (A3 + z * (A5 + z * (A7 + z * (A9 + z * (A11 + z * (A13 + z * (A15 + z * (A17 + z * (A19 + z * A21))))))))));
    else if (a < 6.1227405e-2)
        return t * (1.0 + z * (A3 + z * (A5 + z * (A7 + z * (A9 + z * (A11 + z * (A13 + z * (A15 + z * (A17 + z * (A19 + z * (A21 + z * A23)))))))))));
    else if (a < 7.6353468e-2)
        return t * (1.0 + z * (A3 + z * (A5 + z * (A7 + z * (A9 + z * (A11 + z * (A13 + z * (A15 + z * (A17 + z * (A19 + z * (A21 + z * (A23 + z * A25))))))))))));
    else
    {
        double r = sqrt(-h);
        double y = r * t;

        if (fabs(y) <= 1.0)
            return atanh(y) / r;
        else
            return atanh(1.0 / y) / r;
    }
}

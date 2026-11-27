#include <fenv.h>
#include <iso646.h>
#include <math.h>

#include "cei2.h"
#include "gelbd.h"

void gelbd(double phi, double mc, double *elb, double *eld)
{
    static const double PI = 3.141592653589793238462;

    if (mc >= 0.0)
    {
        if (isfinite(phi))
        {
            if (isinf(mc))
            {
                *elb = 0.0;
                *eld = 0.0;
                return;
            }

            double phix = remainder(phi, PI);
            double n = (phi - phix) / PI;
            phi = fabs(phix);

            double kc = 1.0;
            double sz = 0.0;

            if (mc > 1.0)
            {
                mc = 1.0 / mc;
                kc = sqrt(mc);
                phi = atan(tan(phi) / kc);

                double m = 1.0 - mc;
                double c = cos(phi);
                double x = c * c;
                double d2 = mc + m * x;
                double z = c / sqrt(d2);
                sz = z * sqrt(1.0 - x);
            }

            elbd(phi, 0.5 * PI - phi, mc, elb, eld);
            *elb = copysign(*elb + sz, phix);
            *eld = copysign(*eld - sz, phix);

            if (n == 0.0)
            {
                *elb = *elb * kc;
                *eld = *eld * kc;
                return;
            }
            else
            {
                double celb = ceib(mc);
                double celd = ceid(mc);
                *elb = kc * (*elb + 2.0 * n * celb);
                *eld = kc * (*eld + 2.0 * n * celd);
                return;
            }
        }
        else if (isinf(phi) and isinf(mc))
        {
            feraiseexcept(FE_INVALID);
            *elb = nan("");
            *eld = nan("");
            return;
        }
        else
        {
            *elb = phi;
            *eld = phi;
            return;
        }
    }
    else if (mc < 0.0)
    {
        if (phi == 0.0)
        {
            *elb = 0.0;
            *eld = 0.0;
            return;
        }

        double phix = fabs(phi);

        if (phix < 0.5 * PI)
        {
            double m = 1.0 - mc;
            m = 1.0 / m;
            mc = 1.0 - m;
            double k = sqrt(m);
            phix = asin(sin(phix) / k);

            elbd(phix, 0.5 * PI - phix, mc, elb, eld);
            *elb = copysign(*elb, phi);
            *eld = copysign(*eld, phi);

            *elb = k * (*elb + *eld * mc);
            *eld = *eld * k * m;
            return;
        }
        else if (!isnan(phi))
        {
            feraiseexcept(FE_INVALID);
        }
    }

    *elb = nan("");
    *eld = nan("");
}

void elbd(double phi, double phic, double mc, double *b, double *d)
{
    if (phi < 1.25)
        return elsbd(sin(phi), mc, b, d);

    double m = 1.0 - mc;
    double c = sin(phic);
    double x = c * c;
    double d2 = mc + m * x;

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
        {
            return elcbd(c, mc, b, d);
        }
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
        *b = *b * s;
        *d = *d * s * y;
        return;
    }

    double yy[11];
    double ss[11];
    ss[0] = s;
    int j;

    for (j = 1; j <= 10; ++j)
    {
        y = y / ((1.0 + sqrt(1.0 - y)) * (1.0 + sqrt(1.0 - m * y)));
        yy[j] = y;
        ss[j] = sqrt(y);

        if (y < del)
        {
            serbd(y, m, b, d);
            *b = *b * ss[j];
            *d = *d * ss[j] * y;

            for (int i = 1; i <= j; ++i)
            {
                double sy = ss[j - i] * yy[j - i + 1];
                *b = *b * 2.0 - sy;
                *d = *d * 2.0 + sy;
            }

            return;
        }
    }

    feraiseexcept(FE_INVALID);

    serbd(y, m, b, d);
    *b = *b * ss[j];
    *d = *d * ss[j] * y;

    for (int i = 1; i <= j; ++i)
    {
        double sy = ss[j - i] * yy[j - i + 1];
        *b = *b * 2.0 - sy;
        *d = *d * 2.0 + sy;
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

    double yy[11];
    double ss[11];
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
        {
            s = ss[j];
            elsbd(s, mc, b, dx);

            for (int i = 1; i <= j; ++i)
            {
                double sy = ss[j - i] * yy[j - i + 1];
                *b = *b * 2.0 - sy;
                *dx = *dx * 2.0 + sy;
            }

            return;
        }

        c = sqrt(x);
    }

    feraiseexcept(FE_INVALID);

    s = ss[j];
    elsbd(s, mc, b, dx);

    for (int i = 1; i <= j; ++i)
    {
        double sy = ss[j - i] * yy[j - i + 1];
        *b = *b * 2.0 - sy;
        *dx = *dx * 2.0 + sy;
    }
}

void serbd(double y, double m, double *b, double *d)
{
    static const double F10 = 1.0 / 6.0;
    static const double F20 = 3.0 / 40.0;
    static const double F21 = 2.0 / 40.0;
    static const double F30 = 5.0 / 112.0;
    static const double F31 = 3.0 / 112.0;
    static const double F40 = 35.0 / 1152.0;
    static const double F41 = 20.0 / 1152.0;
    static const double F42 = 18.0 / 1152.0;
    static const double F50 = 63.0 / 2816.0;
    static const double F51 = 35.0 / 2816.0;
    static const double F52 = 30.0 / 2816.0;
    static const double F60 = 231.0 / 13312.0;
    static const double F61 = 126.0 / 13312.0;
    static const double F62 = 105.0 / 13312.0;
    static const double F63 = 100.0 / 13312.0;
    static const double F70 = 429.0 / 30720.0;
    static const double F71 = 231.0 / 30720.0;
    static const double F72 = 189.0 / 30720.0;
    static const double F73 = 175.0 / 30720.0;
    static const double F80 = 6435.0 / 557056.0;
    static const double F81 = 3432.0 / 557056.0;
    static const double F82 = 2722.0 / 557056.0;
    static const double F83 = 2520.0 / 557056.0;
    static const double F84 = 2450.0 / 557056.0;
    static const double F90 = 12155.0 / 1245184.0;
    static const double F91 = 6435.0 / 1245184.0;
    static const double F92 = 5148.0 / 1245184.0;
    static const double F93 = 4620.0 / 1245184.0;
    static const double F94 = 4410.0 / 1245184.0;
    static const double FA0 = 46189.0 / 5505024.0;
    static const double FA1 = 24310.0 / 5505024.0;
    static const double FA2 = 19305.0 / 5505024.0;
    static const double FA3 = 17160.0 / 5505024.0;
    static const double FA4 = 16170.0 / 5505024.0;
    static const double FA5 = 15876.0 / 5505024.0;
    static const double FB0 = 88179.0 / 12058624.0;
    static const double FB1 = 46189.0 / 12058624.0;
    static const double FB2 = 36465.0 / 12058624.0;
    static const double FB3 = 32175.0 / 12058624.0;
    static const double FB4 = 30030.0 / 12058624.0;
    static const double FB5 = 29106.0 / 12058624.0;

    static const double A1 = 3.0 / 5.0;
    static const double A2 = 5.0 / 7.0;
    static const double A3 = 7.0 / 9.0;
    static const double A4 = 9.0 / 11.0;
    static const double A5 = 11.0 / 13.0;
    static const double A6 = 13.0 / 15.0;
    static const double A7 = 15.0 / 17.0;
    static const double A8 = 17.0 / 19.0;
    static const double A9 = 19.0 / 21.0;
    static const double AA = 21.0 / 23.0;
    static const double AB = 23.0 / 25.0;

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

    static const double D0 = 1.0 / 3.0;
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

void rgelbd(float phi, float mc, float *elb, float *eld)
{
    static const float PI = 3.141592653589793238462f;

    if (mc >= 0.0f)
    {
        if (isfinite(phi))
        {
            if (isinf(mc))
            {
                *elb = 0.0f;
                *eld = 0.0f;
                return;
            }

            float phix = remainderf(phi, PI);
            float n = (phi - phix) / PI;
            phi = fabsf(phix);

            float kc = 1.0f;
            float sz = 0.0f;

            if (mc > 1.0f)
            {
                mc = 1.0f / mc;
                kc = sqrtf(mc);
                phi = atanf(tanf(phi) / kc);

                float m = 1.0f - mc;
                float c = cosf(phi);
                float x = c * c;
                float d2 = mc + m * x;
                float z = c / sqrtf(d2);
                sz = z * sqrtf(1.0f - x);
            }

            relbd(phi, 0.5f * PI - phi, mc, elb, eld);
            *elb = copysignf(*elb + sz, phix);
            *eld = copysignf(*eld - sz, phix);

            if (n == 0.0f)
            {
                *elb = *elb * kc;
                *eld = *eld * kc;
                return;
            }
            else
            {
                float celb = sceib(mc);
                float celd = sceid(mc);
                *elb = kc * (*elb + 2.0f * n * celb);
                *eld = kc * (*eld + 2.0f * n * celd);
                return;
            }
        }
        else if (isinf(phi) and isinf(mc))
        {
            feraiseexcept(FE_INVALID);
            *elb = nanf("");
            *eld = nanf("");
            return;
        }
        else
        {
            *elb = phi;
            *eld = phi;
            return;
        }
    }
    else if (mc < 0.0f)
    {
        if (phi == 0.0f)
        {
            *elb = 0.0f;
            *eld = 0.0f;
            return;
        }

        float phix = fabsf(phi);

        if (phix < 0.5 * PI)
        {
            float m = 1.0f - mc;
            m = 1.0f / m;
            mc = 1.0f - m;
            float k = sqrtf(m);
            phix = asinf(sinf(phix) / k);

            relbd(phix, 0.5f * PI - phix, mc, elb, eld);
            *elb = copysignf(*elb, phi);
            *eld = copysignf(*eld, phi);

            *elb = k * (*elb + *eld * mc);
            *eld = *eld * k * m;
            return;
        }
        else if (!isnan(phi))
        {
            feraiseexcept(FE_INVALID);
        }
    }

    *elb = nanf("");
    *eld = nanf("");
}

void relbd(float phi, float phic, float mc, float *b, float *d)
{
    if (phi < 1.25f)
        return relsbd(sinf(phi), mc, b, d);

    float m = 1.0 - mc;
    float c = sinf(phic);
    float x = c * c;
    float d2 = mc + m * x;

    if (x < 0.9f * d2)
    {
        float z = c / sqrtf(d2);
        relsbd(z, mc, b, d);

        float sz = z * sqrtf(1.0f - x);
        *b = sceib(mc) - (*b - sz);
        *d = sceid(mc) - (*d + sz);
    }
    else
    {
        float v = mc * (1.0f - x);

        if (v < x * d2)
        {
            return relcbd(c, mc, b, d);
        }
        else
        {
            float t2 = (1.0f - x) / d2;
            relcbd(sqrtf(mc * t2), mc, b, d);

            float sz = c * sqrtf(t2);
            *b = sceib(mc) - (*b - sz);
            *d = sceid(mc) - (*d + sz);
        }
    }
}

void relsbd(float s, float mc, float *b, float *d)
{
    float m = 1.0f - mc;
    float del = 0.1888f - 0.0378f * m;
    float y = s * s;

    if (y < del)
    {
        rserbd(y, m, b, d);
        *b = *b * s;
        *d = *d * s * y;
        return;
    }

    float yy[11];
    float ss[11];
    ss[0] = s;
    int j;

    for (j = 1; j <= 10; ++j)
    {
        y = y / ((1.0f + sqrtf(1.0f - y)) * (1.0f + sqrtf(1.0f - m * y)));
        yy[j] = y;
        ss[j] = sqrtf(y);

        if (y < del)
        {
            rserbd(y, m, b, d);
            *b = *b * ss[j];
            *d = *d * ss[j] * y;

            for (int i = 1; i <= j; ++i)
            {
                float sy = ss[j - i] * yy[j - i + 1];
                *b = *b * 2.0f - sy;
                *d = *d * 2.0f + sy;
            }

            return;
        }
    }

    feraiseexcept(FE_INVALID);

    rserbd(y, m, b, d);
    *b = *b * ss[j];
    *d = *d * ss[j] * y;

    for (int i = 1; i <= j; ++i)
    {
        float sy = ss[j - i] * yy[j - i + 1];
        *b = *b * 2.0f - sy;
        *d = *d * 2.0f + sy;
    }
}

void relcbd(float c, float mc, float *b, float *dx)
{
    float x = c * c;
    float y = 1.0f - x;
    float s = sqrtf(y);

    if (x > 0.1f)
        return relsbd(s, mc, b, dx);

    float m = 1.0f - mc;

    float yy[11];
    float ss[11];
    ss[0] = s;
    int j;

    for (j = 1; j <= 10; ++j)
    {
        float d = sqrtf(mc + m * x);
        x = (c + d) / (1.0f + d);
        y = 1.0f - x;
        yy[j] = y;
        ss[j] = sqrtf(y);

        if (x > 0.1f)
        {
            s = ss[j];
            relsbd(s, mc, b, dx);

            for (int i = 1; i <= j; ++i)
            {
                float sy = ss[j - i] * yy[j - i + 1];
                *b = *b * 2.0f - sy;
                *dx = *dx * 2.0f + sy;
            }

            return;
        }

        c = sqrtf(x);
    }

    feraiseexcept(FE_INVALID);

    s = ss[j];
    relsbd(s, mc, b, dx);

    for (int i = 1; i <= j; ++i)
    {
        float sy = ss[j - i] * yy[j - i + 1];
        *b = *b * 2.0f - sy;
        *dx = *dx * 2.0f + sy;
    }
}

void rserbd(float y, float m, float *b, float *d)
{
    static const float F10 = 1.0 / 6.0;
    static const float F20 = 3.0 / 40.0;
    static const float F21 = 2.0 / 40.0;
    static const float F30 = 5.0 / 112.0;
    static const float F31 = 3.0 / 112.0;
    static const float F40 = 35.0 / 1152.0;
    static const float F41 = 20.0 / 1152.0;
    static const float F42 = 18.0 / 1152.0;
    static const float F50 = 63.0 / 2816.0;
    static const float F51 = 35.0 / 2816.0;
    static const float F52 = 30.0 / 2816.0;
    static const float F60 = 231.0 / 13312.0;
    static const float F61 = 126.0 / 13312.0;
    static const float F62 = 105.0 / 13312.0;
    static const float F63 = 100.0 / 13312.0;

    static const float A1 = 3.0 / 5.0;
    static const float A2 = 5.0 / 7.0;
    static const float A3 = 7.0 / 9.0;
    static const float A4 = 9.0 / 11.0;
    static const float A5 = 11.0 / 13.0;
    static const float A6 = 13.0 / 15.0;

    float F1 = F10 + m * F10;
    float F2 = F20 + m * (F21 + m * F20);
    float F3 = F30 + m * (F31 + m * (F31 + m * F30));
    float F4 = F40 + m * (F41 + m * (F42 + m * (F41 + m * F40)));
    float F5 = F50 + m * (F51 + m * (F52 + m * (F52 + m * (F51 + m * F50))));
    float F6 = F60 + m * (F61 + m * (F62 + m * (F63 + m * (F62 + m * (F61 + m * F60)))));

    static const float D0 = 1.0 / 3.0;
    float D1 = F1 * A1;
    float D2 = F2 * A2;
    float D3 = F3 * A3;
    float D4 = F4 * A4;
    float D5 = F5 * A5;
    float D6 = F6 * A6;

    *d = D0 + y * (D1 + y * (D2 + y * (D3 + y * (D4 + y * (D5 + y * D6)))));

    float B1 = F1 - D0;
    float B2 = F2 - D1;
    float B3 = F3 - D2;
    float B4 = F4 - D3;
    float B5 = F5 - D4;
    float B6 = F6 - D5;

    *b = 1.0f + y * (B1 + y * (B2 + y * (B3 + y * (B4 + y * (B5 + y * B6)))));
}

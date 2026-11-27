#include <fenv.h>
#include <math.h>

#include "cei2.h"
#include "gelf.h"

double gelf(double phi, double m)
{
    static const double PI = 3.141592653589793238462;

    double mc = 1.0 - m;

    if (mc >= 0.0)
    {
        if (isfinite(phi))
        {
            if (isinf(mc))
                return 0.0;

            double phix = remainder(phi, PI);
            double n = (phi - phix) / PI;
            phi = fabs(phix);

            double kc = 1.0;

            if (mc > 1.0)
            {
                mc = 1.0 / mc;
                kc = sqrt(mc);
                phi = atan(tan(phi) / kc);
            }

            double f = elf(phi, 0.5 * PI - phi, mc);
            f = copysign(f, phix);

            if (n == 0.0)
            {
                return kc * f;
            }
            else
            {
                double k = ceik(mc);
                return kc * (f + 2.0 * n * k);
            }
        }
        else if (isinf(phi) && isinf(mc))
        {
            feraiseexcept(FE_INVALID);
            return nan("");
        }
        else
        {
            return phi;
        }
    }
    else if (mc < 0.0)
    {
        if (phi == 0.0)
            return 0.0;

        double phix = fabs(phi);

        if (phix < 0.5 * PI)
        {
            m = 1.0 / m;
            mc = 1.0 - m;
            double k = sqrt(m);
            phix = asin(sin(phix) / k);

            double f = elf(phix, 0.5 * PI - phix, mc);
            f = copysign(f, phi);

            return k * f;
        }
        else if (!isnan(phi))
        {
            feraiseexcept(FE_INVALID);
        }
    }

    return nan("");
}

double elf(double phi, double phic, double mc)
{
    double m = 1.0 - mc;

    if (phi < 1.249)
        return asn(sin(phi), m);

    double c = sin(phic);
    double yc = c * c;
    double d2 = mc + m * yc;

    if (yc < 0.9 * d2)
        return ceik(mc) - asn(c / sqrt(d2), m);

    double v = mc * (1.0 - yc);

    if (v < yc * d2)
        return acn(c, mc);
    else
        return ceik(mc) - acn(sqrt(v / d2), mc);
}

double acn(double c, double mc)
{
    double m = 1.0 - mc;
    double yc = c * c;

    if (yc > 0.5)
    {
        double y = 1.0 - yc;
        double s = sqrt(y);
        return asn(s, m);
    }

    double f = 1.0;

    for (int i = 1; i <= 20; ++i)
    {
        double d = sqrt(mc + m * yc);
        yc = (c + d) / (1.0 + d);
        f = f * 2.0;

        if (yc > 0.5)
        {
            double y = 1.0 - yc;
            double s = sqrt(y);
            return f * asn(s, m);
        }

        c = sqrt(yc);
    }

    feraiseexcept(FE_INVALID);

    double y = 1.0 - yc;
    double s = sqrt(y);
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
            return f * s * serf(y, m);
        }
    }

    feraiseexcept(FE_INVALID);

    s = sqrt(y);
    return f * s * serf(y, m);
}

double serf(double y, double m)
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

float rgelf(float phi, float m)
{
    static const float PI = 3.141592653589793238462f;

    float mc = 1.0f - m;

    if (mc >= 0.0f)
    {
        if (isfinite(phi))
        {
            if (isinf(mc))
                return 0.0f;

            float phix = remainderf(phi, PI);
            float n = (phi - phix) / PI;
            phi = fabsf(phix);

            float kc = 1.0f;

            if (mc > 1.0f)
            {
                mc = 1.0f / mc;
                kc = sqrtf(mc);
                phi = atanf(tanf(phi) / kc);
            }

            float f = relf(phi, 0.5f * PI - phi, mc);
            f = copysignf(f, phix);

            if (n == 0.0f)
            {
                return kc * f;
            }
            else
            {
                float k = sceik(mc);
                return kc * (f + 2.0f * n * k);
            }
        }
        else if (isinf(phi) && isinf(mc))
        {
            feraiseexcept(FE_INVALID);
            return nanf("");
        }
        else
        {
            return phi;
        }
    }
    else if (mc < 0.0f)
    {
        if (phi == 0.0f)
            return 0.0f;

        float phix = fabsf(phi);

        if (phix < 0.5f * PI)
        {
            m = 1.0f / m;
            mc = 1.0f - m;
            float k = sqrtf(m);
            phix = asinf(sinf(phix) / k);

            float f = relf(phix, 0.5f * PI - phix, mc);
            f = copysignf(f, phi);

            return k * f;
        }
        else if (!isnan(phi))
        {
            feraiseexcept(FE_INVALID);
        }
    }

    return nanf("");
}

float relf(float phi, float phic, float mc)
{
    float m = 1.0f - mc;

    if (phi < 1.249f)
        return rasn(sinf(phi), m);

    float c = sinf(phic);
    float yc = c * c;
    float d2 = mc + m * yc;

    if (yc < 0.9f * d2)
        return sceik(mc) - rasn(c / sqrtf(d2), m);

    float v = mc * (1.0f - yc);

    if (v < yc * d2)
        return racn(c, mc);
    else
        return sceik(mc) - racn(sqrtf(v / d2), mc);
}

float racn(float c, float mc)
{
    float m = 1.0f - mc;
    float yc = c * c;

    if (yc > 0.5)
    {
        float y = 1.0f - yc;
        float s = sqrtf(y);
        return rasn(s, m);
    }

    float f = 1.0f;

    for (int i = 1; i <= 10; ++i)
    {
        float d = sqrtf(mc + m * yc);
        yc = (c + d) / (1.0f + d);
        f = f * 2.0f;

        if (yc > 0.5f)
        {
            float y = 1.0f - yc;
            float s = sqrtf(y);
            return f * rasn(s, m);
        }

        c = sqrtf(yc);
    }

    feraiseexcept(FE_INVALID);

    float y = 1.0f - yc;
    float s = sqrtf(y);
    return f * rasn(s, m);
}

float rasn(float s, float m)
{
    float del = 0.1888f - 0.0378f * m;
    float y = s * s;

    if (y < del)
        return s * rserf(y, m);

    float f = 1.0f;

    for (int j = 1; j <= 10; ++j)
    {
        y = y / ((1.0f + sqrtf(1.0f - y)) * (1.0f + sqrtf(1.0f - m * y)));
        f = f * 2.0f;

        if (y < del)
        {
            s = sqrtf(y);
            return f * s * rserf(y, m);
        }
    }

    feraiseexcept(FE_INVALID);

    s = sqrtf(y);
    return f * s * serf(y, m);
}

float rserf(float y, float m)
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

    float F1 = F10 + m * F10;
    float F2 = F20 + m * (F21 + m * F20);
    float F3 = F30 + m * (F31 + m * (F31 + m * F30));
    float F4 = F40 + m * (F41 + m * (F42 + m * (F41 + m * F40)));
    float F5 = F50 + m * (F51 + m * (F52 + m * (F52 + m * (F51 + m * F50))));
    float F6 = F60 + m * (F61 + m * (F62 + m * (F63 + m * (F62 + m * (F61 + m * F60)))));

    return 1.0f + y * (F1 + y * (F2 + y * (F3 + y * (F4 + y * (F5 + y * F6)))));
}

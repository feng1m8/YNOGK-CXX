#include <fenv.h>
#include <iso646.h>
#include <math.h>

#include "cei2.h"
#include "gscd.h"

void gscd(double u, double mc, double *s, double *c, double *d)
{
    const double PI = 3.141592653589793238462;

    if (u <= 0.25 * PI)
    {
        scd2(u, mc, s, c, d);
        return;
    }
    else
    {
        feraiseexcept(FE_INVALID);

        *s = nan("");
        *c = nan("");
        *d = nan("");
        return;
    }
}

void scd2(double u, double mc, double *s, double *c, double *d)
{
    static const double B10 = 1.0 / 24.0;
    static const double B11 = 1.0 / 6.0;
    static const double B20 = 1.0 / 720.0;
    static const double B21 = 11.0 / 180.0;
    static const double B22 = 1.0 / 45.0;

    double m = 1.0 - mc;
    double uT = 5.217e-3 - m * 2.143e-3;
    double u0 = u;

    int n;
    for (n = 0; n <= 20; ++n)
    {
        if (u0 < uT)
            break;

        u0 = u0 * 0.5;
    }

    if (n > 20)
        feraiseexcept(FE_INVALID);

    double uA = 1.76269 + mc * 1.16357;
    double v = u0 * u0;
    double a = 1.0;
    double b = v * (0.5 - v * (B10 + m * B11 - v * (B20 + m * (B21 + m * B22))));

    if (u < uA)
    {
        for (int j = 1; j <= n; ++j)
        {
            double y = b * (2.0 * a - b);
            double z = a * a;
            double my = m * y;
            b = 2.0 * y * (z - my);
            a = z * z - my * y;
        }
    }
    else
    {
        for (int j = 1; j <= n; ++j)
        {
            double y = b * (2.0 * a - b);
            double z = a * a;
            double my = m * y;

            if (z < 2.0 * my)
            {
                *c = a - b;
                double mc2 = 2.0 * mc;
                double m2 = 2.0 * m;

                for (int i = j; i <= n; ++i)
                {
                    double x = *c * *c;
                    z = a * a;
                    double w = m * x * x - mc * z * z;
                    double xz = x * z;
                    *c = mc2 * xz + w;
                    a = m2 * xz - w;
                }

                *c = *c / a;
                double x = *c * *c;
                *s = sqrt(1.0 - x);
                *d = sqrt(mc + m * x);
                return;
            }

            b = 2.0 * y * (z - my);
            a = z * z - my * y;
        }
    }

    b = b / a;
    double y = b * (2.0 - b);
    *c = 1.0 - b;
    *s = sqrt(y);
    *d = sqrt(1.0 - m * y);
}

void rscd2(float u, float mc, float *s, float *c, float *d)
{
    float m = 1.0f - mc;
    float uT = 9.207e-4f - m * 4.604e-4f;
    float u0 = u;

    int n;
    for (n = 0; n <= 20; ++n)
    {
        if (u0 < uT)
            break;

        u0 = u0 * 0.5f;
    }

    if (n > 20)
        feraiseexcept(FE_INVALID);

    float uA = 1.76269f + mc * 1.16357f;
    float v = u0 * u0;
    float a = 1.0f;
    float b = v * 0.5f;

    if (u < uA)
    {
        for (int j = 1; j <= n; ++j)
        {
            float y = b * (2.0f * a - b);
            float z = a * a;
            float my = m * y;
            b = 2.0f * y * (z - my);
            a = z * z - my * y;
        }
    }
    else
    {
        for (int j = 1; j <= n; ++j)
        {
            float y = b * (2.0f * a - b);
            float z = a * a;
            float my = m * y;

            if (z < 2.0f * my)
            {
                *c = a - b;
                float mc2 = 2.0f * mc;
                float m2 = 2.0f * m;

                for (int i = j; i <= n; ++i)
                {
                    float x = *c * *c;
                    z = a * a;
                    float w = m * x * x - mc * z * z;
                    float xz = x * z;
                    *c = mc2 * xz + w;
                    a = m2 * xz - w;
                }

                *c = *c / a;
                float x = *c * *c;
                *s = sqrtf(1.0f - x);
                *d = sqrtf(mc + m * x);
                return;
            }

            b = 2.0f * y * (z - my);
            a = z * z - my * y;
        }
    }

    b = b / a;
    float y = b * (2.0f - b);
    *c = 1.0f - b;
    *s = sqrtf(y);
    *d = sqrtf(1.0f - m * y);
}

#include <algorithm>
#include <complex>
#include <cstring>
#include <functional>
#include <iterator>

extern "C"
{
#include "solve2345/solve2345.h"
}

extern "C" void sort(double a1, double a2, double a3, double *s)
{
    double temp;
    int i, j;

    s[0] = a1;
    s[1] = a2;
    s[2] = a3;

    for (i = 0; i < 3; i++)
    {
        for (j = i + 1; j < 3; j++)
        {
            if (s[i] < s[j])
            {
                temp = s[i];
                s[i] = s[j];
                s[j] = temp;
            }
        }
    }
}

extern "C" int root3(double b, double c, double d, double _Complex *roots, unsigned int *reals)
{
    double zeror[3], zeroi[3];

    if (b == 0.0)
        cubicSolveDepressed(c, d, zeror, zeroi);
    else
        cubicSolveMonic(b, c, d, zeror, zeroi);

    if (zeroi[0] == 0.0)
    {
        std::sort(std::begin(zeror), std::end(zeror), std::greater<double>());
        roots[1] = zeror[0];
        roots[2] = zeror[1];
        roots[3] = zeror[2];

        *reals = 3;
    }
    else
    {
        std::complex<double> roots_cxx[3]{
            {zeror[2], zeroi[2]},
            {zeror[0], zeroi[0]},
            {zeror[1], zeroi[1]},
        };

        static_assert(sizeof(std::complex<double>) == sizeof(double _Complex));
        std::memcpy(roots + 1, roots_cxx, sizeof(roots_cxx));

        *reals = 1;
    }

    return 0;
}

extern "C" int root4(double b, double c, double d, double e, double _Complex *roots, unsigned int *reals)
{
    double zeror[4], zeroi[4];

    if (b == 0.0)
        quarticSolveDepressed(c, d, e, zeror, zeroi);
    else
        quarticSolveMonic(b, c, d, e, zeror, zeroi);

    std::complex<double> roots_cxx[4];
    *reals = 0;
    for (int i = 0; i < 4; ++i)
    {
        roots_cxx[i] = {zeror[i], zeroi[i]};

        if (zeroi[i] == 0.0)
            ++(*reals);
    }

    std::sort(
        std::begin(roots_cxx), std::end(roots_cxx),
        [](const auto &x, const auto &y)
        {
            return x.real() < y.real();
        });

    static_assert(sizeof(std::complex<double>) == sizeof(double _Complex));
    std::memcpy(roots + 1, roots_cxx, sizeof(roots_cxx));

    return 0;
}

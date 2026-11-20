#include <cmath>

#include "elbd.h"
#include "elbdj2.h"
#include "elf.h"

extern "C"
{
#include "root34.h"

    double rc(double x, double y)
    {
        if (y == 0.0)
            return infinity;

        double t = 1.0 / std::sqrt(x);
        double h = y - x;

        if (y < 0.0)
        {
            double ri = 1.0 / std::sqrt(-h);
            return uatan(ri, -x) * ri / t;
        }

        return uatan(t, h);
    }

    double rf(double x, double y, double z)
    {
        if (x < 0.0)
            x = 0.0;
        if (y < 0.0)
            y = 0.0;
        if (z < 0.0)
            z = 0.0;

        double s[3];
        sort(x, y, z, s);

        double phi = std::acos(std::sqrt(s[2] / s[0]));
        double mc = (s[1] - s[2]) / (s[0] - s[2]);
        return elf(phi, 0.5 * pi - phi, mc) / std::sqrt(s[0] - s[2]);
    }

    double rd(double x, double y, double z)
    {
        if (x < 0.0)
            x = 0.0;
        if (y < 0.0)
            y = 0.0;

        double s[3];
        sort(x, y, z, s);

        double temp = std::sqrt(s[0] - s[2]);

        if (z == s[1])
        {
            if (s[0] == s[1])
            {
                double phi = std::acos(std::sqrt(s[2] / s[0]));
                double t1 = temp * std::sqrt(s[2]) / s[0];
                return 1.5 * (phi - t1) / (temp * (s[0] - s[2]));
            }
            if (s[1] == s[2])
            {
                double t1 = std::sqrt(s[0]) / (temp * s[2]);
                return 1.5 * (t1 - rc(s[0], s[2])) / temp;
            }
        }

        double phi = std::acos(std::sqrt(s[2] / s[0]));
        double mc = (s[1] - s[2]) / (s[0] - s[2]);
        double elb, eld;
        elbd(phi, 0.5 * pi - phi, mc, &elb, &eld);

        if (z == s[0])
            return 3.0 * eld / (temp * (s[0] - s[2]));
        else if (z == s[1])
        {
            double t1 = temp * std::sqrt(s[2] / (s[0] * s[1]));
            return 3.0 * (elb - t1) / (temp * (s[1] - s[2]));
        }
        else
        {
            double el2 = elb + mc * eld;
            double t1 = temp * std::sqrt(s[1] / (s[0] * s[2]));
            return 3.0 * (t1 - el2) / (temp * (s[1] - s[2]));
        }
    }

    double rj(double x, double y, double z, double p)
    {
        if (p == x or p == y or p == z)
            return rd(x, y, z);

        if (x < 0.0)
            x = 0.0;
        if (y < 0.0)
            y = 0.0;
        if (z < 0.0)
            z = 0.0;

        double s[3];
        sort(x, y, z, s);

        double phi = std::acos(std::sqrt(s[2] / s[0]));
        double m = (s[0] - s[1]) / (s[0] - s[2]);
        double mc = 1.0 - m;

        double elb, eld, elj;

        double temp = std::sqrt(s[0] - s[2]);
        if (p <= 0.0)
        {
            double a = 1.0 / (s[1] - p);
            double b = a * (s[0] - s[1]) * (s[1] - s[2]);
            double pt = s[1] + b;
            double rho = s[0] * s[2] / s[1];
            double tau = p * pt / s[1];

            double n = (s[0] - pt) / (s[0] - s[2]);

            elbdj2(phi, 0.5 * pi - phi, n, mc, &elb, &eld, &elj);
            return 3.0 * a * (rc(rho, tau) + ((m - n) * elj - elb - eld) / temp);
        }

        double n = (s[0] - p) / (s[0] - s[2]);
        if (n > 1.0)
        {
            double n1 = m / n;
            double t1 = temp * std::sqrt(s[0] / (s[1] * s[2]));
            double h1 = (1.0 - n) * (1.0 - n1);
            elbdj2(phi, 0.5 * pi - phi, n1, mc, &elb, &eld, &elj);
            elj = (-elb - eld + uatan(t1, h1) - n1 * elj) / n;
        }
        else if (n < 0.0)
        {
            double nc = 1.0 - n;
            double n1c = mc / nc;
            double n1 = 1.0 - n1c;
            double t1 = temp * std::sqrt(s[2] / (s[0] * s[1]));
            double h1 = -n * n1;
            elbdj2(phi, 0.5 * pi - phi, n1, mc, &elb, &eld, &elj);
            elj = (elb + eld - uatan(t1, h1) - n1c * elj) / nc;
        }
        else
            elbdj2(phi, 0.5 * pi - phi, n, mc, &elb, &eld, &elj);

        return 3.0 * elj / (temp * (s[0] - s[2]));
    }
}

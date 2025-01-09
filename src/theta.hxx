#ifndef YNOGK_CXX_THETA_HXX
#define YNOGK_CXX_THETA_HXX

#include <variant>

extern "C"
{
#include "ellCarlsons.h"
#include "ynogk_theta_part.h"
}

namespace ynogk
{
    class ThetaKerr
    {
        friend class Theta;

    private:
#define static
#include "ynogk_theta_part.c"
#undef static

        int Get_Integrations_of_Theta_part(ptcl *pt, double p, double *phit, double *mucos, double *sign_pth, int *tm1, int *tm2);
    };

    class ThetaSchwarzschild
    {
        friend class Theta;

    private:
#define static
#include "ynogk_0aspin.c"
#undef static

        int Integration_Theta_part_Settings(ptcl *pt)
        {
            this->phit_Schwarzschild_Settings(pt);
            return 0;
        }

        int Get_Integrations_of_Theta_part(ptcl *pt, double p, double *phit, double *mucos, double *sign_pth, int *tm1, int *tm2)
        {
            return this->Get_phit_Schwarzschild(pt, p, phit, mucos, sign_pth, tm1, tm2);
        }

        void Get_Integrals_For_Theta_Part(ptcl *pt, double p, double *phit, double *timet, double *mucos, double *sign_pth, int *tm1, int *tm2)
        {
            *timet = 0.0;
            this->Get_phit_Integrals_Schwarzschild(pt, p, phit, mucos, sign_pth, tm1, tm2);
        }
    };

    class Theta
    {
        friend class Particle;

    private:
        std::variant<ThetaKerr, ThetaSchwarzschild> theta;

        Theta(double a_spin)
        {
            if (a_spin == 0.0)
                this->theta.emplace<ThetaSchwarzschild>();
        }

        int Integration_Theta_part_Settings(ptcl *pt)
        {
            return std::visit(
                [&](auto &self)
                {
                    return self.Integration_Theta_part_Settings(pt);
                },
                this->theta);
        }

        int Get_Integrations_of_Theta_part(ptcl *pt, double p, double *phit, double *mucos, double *sign_pth, int *tm1, int *tm2)
        {
            return std::visit(
                [&](auto &self)
                {
                    return self.Get_Integrations_of_Theta_part(pt, p, phit, mucos, sign_pth, tm1, tm2);
                },
                this->theta);
        }

        void Get_Integrals_For_Theta_Part(ptcl *pt, double p, double *phit, double *timet, double *mucos, double *sign_pth, int *tm1, int *tm2)
        {
            std::visit(
                [&](auto &self)
                {
                    return self.Get_Integrals_For_Theta_Part(pt, p, phit, timet, mucos, sign_pth, tm1, tm2);
                },
                this->theta);
        }
    };
}

#endif

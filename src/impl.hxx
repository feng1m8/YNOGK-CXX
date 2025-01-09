#ifndef YNOGK_CXX_IMPL_HXX
#define YNOGK_CXX_IMPL_HXX

#include "ynogk_cxx/particle.hxx"

#include "radius.hxx"
#include "theta.hxx"

extern "C"
{
#include "ynogkBL.h"
}

namespace ynogk
{
    class Particle::Impl
    {
        friend class Particle;

    private:
#include "ynogkBL.c"

        Radius radius;
        Theta theta;

        Impl(double a_spin) : theta(a_spin){};

        double phi(ptcl *p, double pm);

        void Get_Integrals_For_R_Part(ptcl *pt, double p, double *r_coord, double *sign_pr, double *affr, double *timer, double *phir)
        {
            this->radius.Get_Integrals_For_R_Part(pt, p, r_coord, sign_pr, affr, timer, phir);
        }

        void Get_Integrals_For_Theta_Part(ptcl *pt, double p, double *phit, double *timet, double *mucos, double *sign_pth, int *tm1, int *tm2)
        {
            this->theta.Get_Integrals_For_Theta_Part(pt, p, phit, timet, mucos, sign_pth, tm1, tm2);
        }
    };
}

#endif

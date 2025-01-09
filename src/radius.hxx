#ifndef YNOGK_CXX_RADIUS_HXX
#define YNOGK_CXX_RADIUS_HXX

#include "coordinates.hxx"

extern "C"
{
#include "ynogkini.h"
}

namespace ynogk
{
    class Radius
    {
        friend class Particle;

    private:
#define static
#include "ynogkini.c"
#undef static

        Coordinates coordinates;

        void Integral_r_part(ptcl *p, double pm, double *phir);
        void Get_Results_of_Integrations_For_R_Part(ptcl *p, double pm, double *r_coord, double *sign_pr, double *phir);

        double radius_settings(ptcl *p)
        {
            return this->coordinates.radius_settings(p);
        }

        void Get_data(out_data2 *tmp)
        {
            this->coordinates.Get_data(tmp);
        }
    };
}

#endif

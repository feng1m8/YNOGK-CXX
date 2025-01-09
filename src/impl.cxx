#include "impl.hxx"

namespace ynogk
{
    double Particle::Impl::phi(ptcl *p, double pm)
    {
        double phi_r, phi_t, Rab;
        int tm1, tm2;

        this->radius.Set_Initializations_For_Integrations_of_R_Part(p);
        this->radius.Get_Results_of_Integrations_For_R_Part(p, pm, &p->r_p, &p->sign_pr_p, &phi_r);

        this->theta.Integration_Theta_part_Settings(p);
        this->theta.Get_Integrations_of_Theta_part(p, pm, &phi_t, &p->mu_p, &p->sign_pth_p, &tm1, &tm2);

        p->sin_p = sqrt(one - p->mu_p * p->mu_p);

        if (fabs(p->muobs) < one)
        {
            p->phi_p = -(phi_r + phi_t);

            if (p->f1234[3] == zero)
            {
                if (p->mu_tp1 == one)
                    p->phi_p = p->phi_p + tm1 * pi;

                if (p->mu_tp2 == -one)
                    p->phi_p = p->phi_p + tm2 * pi;
            }
        }
        else
        {
            p->phi_p = -(phi_t + phi_r);

            if (p->mu_tp1 == one)
                p->phi_p = p->phi_p + tm1 * pi;

            if (p->mu_tp2 == -one)
                p->phi_p = p->phi_p + tm2 * pi;

            Rab = sqrt(sq(p->f1234[3]) + sq(p->f1234[2]));

            if (Rab > zero)
            {
                if ((p->f1234[3] >= zero) && (p->f1234[2] > zero))
                    p->phi_p = p->muobs * p->phi_p + asin(p->f1234[2] / Rab);

                if ((p->f1234[3] < zero) && (p->f1234[2] >= zero))
                    p->phi_p = p->muobs * p->phi_p + pi - asin(p->f1234[2] / Rab);

                if ((p->f1234[3] <= zero) && (p->f1234[2] < zero))
                    p->phi_p = p->muobs * p->phi_p + pi - asin(p->f1234[2] / Rab);

                if ((p->f1234[3] > zero) && (p->f1234[2] <= zero))
                    p->phi_p = p->muobs * p->phi_p + twopi + asin(p->f1234[2] / Rab);
            }
            else
                p->phi_p = zero;
        }

        return p->phi_p;
    }
}

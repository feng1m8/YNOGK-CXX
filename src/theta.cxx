#include "theta.hxx"

namespace ynogk
{
    int ThetaKerr::Get_Integrations_of_Theta_part(ptcl *pt, double p, double *phit, double *mucos, double *sign_pth, int *tm1, int *tm2)
    {
        double pp_phi, p1_phi, p2_phi, tmu, u;

        *tm1 = 0;
        *tm2 = 0;

        if (pt->f1234[3] == zero && pt->f1234[2] == zero && fabs(muobs) == one)
        {
            *sign_pth = zero;
            *mucos = sign(one, muobs);
            *phit = zero;
            return 0;
        }

        if (muobs == zero && (fabs(lambda) < fabs(a_spin)) && q == zero)
        {
            *sign_pth = zero;
            *phit = zero;
            *mucos = zero;
            return 0;
        }

        if (fabs(mu_tp1 - zero) < 1.e-11)
        {
            *sign_pth = zero;
            *phit = zero;
            *mucos = zero;
            *tm1 = 0;
            *tm2 = 0;
            return 0;
        }

        if (a_spin == zero)
            return 0;

        tmu = weierstrassP(p + PI01, g2, g3, dd, del);
        *mucos = mu_tp1 + b0 / (four * tmu - b1);

        u = p + PI01;
        if (u <= zero)
            *sign_pth = -one;
        else
        {
            u = fmod(u, period_wp);
            if (u <= half_period_wp)
                *sign_pth = one;
            else
                *sign_pth = -one;
        }

        if (lambda != zero)
        {
            index_p4[1] = 0;
            cases_int = 1;
            weierstrass_int_J3(tobs, tmu, dd, del, a4, b4, index_p4, rff_p, integ4, cases_int);
            pp = integ4[1];
        }

        Get_mu_t1_t2_New(p, pt->f1234[2], mobseqmtp, p_mt1_mt2, muobs, mu_tp1, mu_tp2, &t1, &t2);
        *tm1 = t1;
        *tm2 = t2;

        index_p4[1] = -1;
        index_p4[2] = -2;
        index_p4[3] = 0;
        index_p4[4] = -4;

        if (lambda != zero)
        {
            cases_int = 2;
            weierstrass_int_J3(tobs, tmu, dd, del, -tplus, b4, index_p4, fabs(pp), integ4, cases_int);
            weierstrass_int_J3(tobs, tmu, dd, del, -tminus, b4, index_p4, fabs(pp), integ14, cases_int);
            pp_phi = lambda * (pp / (one - mu_tp12) + integ4[2] * Wmup - integ14[2] * Wmum);
        }
        else
            pp_phi = zero;

        if (t1 == 0)
            p1_phi = zero;
        else
        {
            if (lambda != zero)
            {
                cases_int = 2;
                weierstrass_int_J3(tobs, infinity, dd, del, -tplus, b4, index_p4, PI0, integ4, cases_int);
                weierstrass_int_J3(tobs, infinity, dd, del, -tminus, b4, index_p4, PI0, integ14, cases_int);
                PI1_phi = lambda * (PI0 / (one - mu_tp12) + integ4[2] * Wmup - integ14[2] * Wmum);

                p1_phi = PI1_phi - pp_phi;
            }
            else
                p1_phi = zero;
        }

        if (t2 == 0)
            p2_phi = zero;
        else
        {
            if (lambda != zero)
            {
                cases_int = 2;
                weierstrass_int_J3(tp2, tobs, dd, del, -tplus, b4, index_p4, PI2_p, integ4, cases_int);
                weierstrass_int_J3(tp2, tobs, dd, del, -tminus, b4, index_p4, PI2_p, integ14, cases_int);
                PI2_phi = lambda * (PI2_p / (one - mu_tp12) + integ4[2] * Wmup - integ14[2] * Wmum);

                p2_phi = PI2_phi + pp_phi;
            }
            else
                p2_phi = zero;
        }

        if (mobseqmtp)
        {
            if (muobs == mu_tp1)
                *phit = -pp_phi + two * (t1 * p1_phi + t2 * p2_phi);
            else
                *phit = pp_phi + two * (t1 * p1_phi + t2 * p2_phi);
        }
        else
        {
            if (pt->f1234[2] < zero)
                *phit = pp_phi + two * (t1 * p1_phi + t2 * p2_phi);

            if (pt->f1234[2] > zero)
                *phit = -pp_phi + two * (t1 * p1_phi + t2 * p2_phi);
        }

        return 0;
    }
}

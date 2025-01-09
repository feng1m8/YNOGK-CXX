#include "radius.hxx"

namespace ynogk
{
    void Radius::Integral_r_part(ptcl *p, double pm, double *phir)
    {
        double pp_phi, p1_phi, p2_phi;

        index_p4[1] = -1;
        index_p4[2] = -2;
        index_p4[3] = 0;
        index_p4[4] = -4;

        b4 = 1.0;

        if (p->a_spin != zero)
        {
            cases_int = 2;
            weierstrass_int_J3(tobs, tp, dd, del, -E_add, b4, index_p4, fabs(pp), integ4, cases_int);
            weierstrass_int_J3(tobs, tp, dd, del, -E_m, b4, index_p4, fabs(pp), integ14, cases_int);
            pp_phi = pp * p->a_spin * (B_add / (r_tp1 - r_add) - B_m / (r_tp1 - r_m)) - p->a_spin * B_add * D_add * integ4[2] + p->a_spin * B_m * D_m * integ14[2];
        }
        else
            pp_phi = zero;

        if (p->muobs == zero && p->f1234[2] == zero)
            pp_phi = pp_phi + pp * p->lambda;

        if (t1 == 0)
            p1_phi = zero;
        else
        {
            if (p->a_spin != zero)
            {
                cases_int = 2;
                weierstrass_int_J3(tobs, infinity, dd, del, -E_add, b4, index_p4, PI0, integ4, cases_int);
                weierstrass_int_J3(tobs, infinity, dd, del, -E_m, b4, index_p4, PI0, integ14, cases_int);
                PI1_phi = PI0 * p->a_spin * (B_add / (r_tp1 - r_add) - B_m / (r_tp1 - r_m)) - p->a_spin * B_add * D_add * integ4[2] + p->a_spin * B_m * D_m * integ14[2];
            }
            else
                PI1_phi = zero;

            if (p->muobs == zero && p->f1234[2] == zero)
                PI1_phi += PI0 * lambda;

            p1_phi = PI1_phi - pp_phi;
        }

        if (t2 == zero)
            p2_phi = zero;
        else
        {
            if (p->a_spin != zero)
            {
                cases_int = 2;
                weierstrass_int_J3(tp2, tobs, dd, del, -E_add, b4, index_p4, PI2_p, integ4, cases_int);
                weierstrass_int_J3(tp2, tobs, dd, del, -E_m, b4, index_p4, PI2_p, integ14, cases_int);
                PI2_phi = PI2_p * p->a_spin * (B_add / (r_tp1 - r_add) - B_m / (r_tp1 - r_m)) - p->a_spin * B_add * D_add * integ4[2] + p->a_spin * B_m * D_m * integ14[2];
            }
            else
                PI2_phi = zero;

            if (p->muobs == zero && p->f1234[2] == zero)
                PI2_phi = PI2_phi + PI2_p * p->lambda;

            p2_phi = PI2_phi + pp_phi;
        }

        if (p->f1234[1] != zero)
            *phir = sign(one, -p->f1234[1]) * pp_phi + two * (t1 * p1_phi + t2 * p2_phi);
        else
        {
            if (robs == r_tp1)
                *phir = -pp_phi + two * (t1 * p1_phi + t2 * p2_phi);
            else
                *phir = pp_phi + two * (t1 * p1_phi + t2 * p2_phi);
        }
    }

    void Radius::Get_Results_of_Integrations_For_R_Part(ptcl *p, double pm, double *r_coord, double *sign_pr, double *phir)
    {
        if (p->r_reals != 0)
        {
            if (r_cases == 1)
            {
                if (p->f1234[1] >= zero)
                {
                    if (pm < PI0_obs_inf)
                    {
                        tp = weierstrassP(pm + PI0, g2, g3, dd, del);
                        *r_coord = r_tp1 + b0 / (four * tp - b1);
                        pp = -pm;
                    }
                    else
                    {
                        tp = tinf;
                        *r_coord = infinity;
                        pp = -PI0_obs_inf;
                    }

                    t1 = 0;
                    t2 = 0;
                    *sign_pr = one;
                }
                else
                {
                    if (!indrhorizon)
                    {
                        t2 = 0;

                        if (pm <= PI0)
                        {
                            t1 = 0;
                            pp = pm;
                            tp = weierstrassP(pm - PI0, g2, g3, dd, del);
                            *r_coord = r_tp1 + b0 / (four * tp - b1);
                            *sign_pr = -one;
                        }
                        else
                        {
                            t1 = 1;
                            PI1_p = PI0;

                            if (pm < PI0_total)
                            {
                                tp = weierstrassP(pm - PI0, g2, g3, dd, del);
                                *r_coord = r_tp1 + b0 / (four * tp - b1);
                                pp = two * PI0 - pm;
                            }
                            else
                            {
                                tp = tinf;
                                *r_coord = infinity;
                                pp = -PI0_total + two * PI0;
                            }

                            *sign_pr = one;
                        }
                    }
                    else
                    {
                        if (pm < PI0_obs_hori)
                        {
                            tp = weierstrassP(pm - PI0, g2, g3, dd, del);
                            *r_coord = r_tp1 + b0 / (four * tp - b1);
                            pp = pm;
                        }
                        else
                        {
                            tp = thorizon;
                            *r_coord = rhorizon;
                            pp = PI0_obs_hori;
                        }

                        t1 = 0;
                        t2 = 0;
                        *sign_pr = -one;
                    }
                }
            }
            else if (r_cases == 2)
            {
                if (!indrhorizon)
                {
                    tp = weierstrassP(pm + PI01, g2, g3, dd, del);
                    *r_coord = r_tp1 + b0 / (four * tp - b1);
                    weierstrass_int_J3(tobs, tp, dd, del, zero, one, index_p4, rff_p, integ4, cases_int);
                    pp = integ4[1];
                    Get_t1_t2(pm, p->f1234[1], robs_eq_rtp, p_tp1_tp2, p->robs, r_tp1, r_tp2, &t1, &t2);
                }
                else
                {
                    if (p->f1234[1] <= zero)
                    {
                        if (pm < PI0_obs_hori)
                        {
                            tp = weierstrassP(pm - PI0, g2, g3, dd, del);
                            *r_coord = r_tp1 + b0 / (four * tp - b1);
                            pp = pm;
                        }
                        else
                        {
                            tp = thorizon;
                            *r_coord = rhorizon;
                            pp = PI0_obs_hori;
                        }

                        t1 = 0;
                        t2 = 0;
                        *sign_pr = -one;
                    }
                    else
                    {
                        if (pm <= PI0_obs_tp2)
                        {
                            t1 = 0;
                            t2 = 0;
                            pp = -pm;
                            tp = weierstrassP(pm + PI0, g2, g3, dd, del);
                            *r_coord = r_tp1 + b0 / (four * tp - b1);
                            *sign_pr = one;
                        }
                        else
                        {
                            t1 = 0;
                            t2 = 1;
                            if (pm < PI0_total)
                            {
                                tp = weierstrassP(pm + PI0, g2, g3, dd, del);
                                *r_coord = r_tp1 + b0 / (four * tp - b1);
                                pp = pm - two * PI0_obs_tp2;
                            }
                            else
                            {
                                tp = thorizon;
                                *r_coord = rhorizon;
                                pp = PI0_total - two * PI0_obs_tp2;
                            }
                            *sign_pr = -one;
                        }
                    }
                }
            }
            Integral_r_part(p, pm, phir);
        }
        else
        {
            if (u != zero)
            {
                if (pm < PI0)
                    sncndn(pm * w * sqrt_L1 + sign(one, p->f1234[1]) * pinf * w * sqrt_L1, one - m2, &sn, &cn, &dn);

                if (p->f1234[1] < zero)
                {
                    if (pm < PI0)
                    {
                        y = u + (-two * u + w * (L1 - L2) * sn * fabs(cn)) / ((L1 - L2) * sq(sn) - (L1 - one));
                        *r_coord = y;
                        pp = pm;
                    }
                    else
                    {
                        y = rhorizon;
                        *r_coord = y;
                        pp = PI0;
                    }

                    x = robs;
                    *sign_pr = -one;
                }
                else
                {
                    if (pm < PI0)
                    {
                        x = u + (-two * u - w * (L1 - L2) * sn * fabs(cn)) / ((L1 - L2) * sq(sn) - (L1 - one));
                        *r_coord = x;
                        pp = pm;
                    }
                    else
                    {
                        x = infinity;
                        *r_coord = x;
                        pp = PI0;
                    }

                    y = robs;
                    *sign_pr = one;
                }

                cases_int = 2;
                carlson_doublecomplex5(y, x, f1, g1, h1, f2, g2, h2, -r_add, b5, index_p5, fabs(pp), integ5, cases_int);
                carlson_doublecomplex5(y, x, f1, g1, h1, f2, g2, h2, -r_m, b5, index_p5, fabs(pp), integ15, cases_int);
                *phir = a_spin * (B_add * integ5[2] - B_m * integ15[2]);

                if (muobs == zero && p->f1234[2] == zero)
                    *phir = *phir + pp * lambda;
            }
            else
            {
                if (p->f1234[1] < zero)
                {
                    PI0 = (atan(robs / w) - atan(rhorizon / w)) / w;

                    if (pm < PI0)
                    {
                        y = w * tan(atan(robs / w) - pm * w);
                        *r_coord = y;
                    }
                    else
                    {
                        y = rhorizon;
                        *r_coord = y;
                    }

                    x = robs;
                    *sign_pr = -one;
                }
                else
                {
                    PI0 = (halfpi - atan(robs / w)) / w;

                    if (pm < PI0)
                    {
                        x = w * tan(atan(robs / w) + pm * w);
                        *r_coord = x;
                    }
                    else
                    {
                        x = infinity;
                        *r_coord = x;
                    }

                    y = robs;
                    *sign_pr = one;
                }

                if (a_spin != zero)
                {
                    *phir = (-B_add * r_add / w / (r_add2 + w2) + B_m * r_m / w / (r_m2 + w2)) * (atan(x / w) - atan(y / w)) +
                            log(fabs(x - r_add) / sqrt(x * x + w2)) * B_add / (r_add * two + w2) -
                            log(fabs(y - r_add) / sqrt(y * y + w2)) * B_add / (r_add * two + w2) -
                            log(fabs(x - r_m) / sqrt(x * x + w2)) * B_m / (r_m * two + w2) +
                            log(fabs(y - r_m) / sqrt(y * y + w2)) * B_m / (r_m * two + w2);
                    *phir = *phir * a_spin;
                }
                else
                    *phir = zero;

                if (muobs == zero && p->f1234[2] == zero)
                    *phir = *phir + lambda * (atan(x / w) - atan(y / w)) / w;
            }
        }
    }
}

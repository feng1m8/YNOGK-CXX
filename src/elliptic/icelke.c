#include <fenv.h>
#include <math.h>

#include "rootfind/zeros.h"

#include "cei2.h"

static double icelk_func(double mc, void *K)
{
    double k = *(double *)K;

    if (mc == 0.0)
        return -k;
    else
        return sqrt(mc) * ceik(mc) - k;
}

double icelk(double K)
{
    static const double KA = 1.5707963267948966;
    static const double KB = 1.7539969906494259;
    static const double KC = 1.9495677498060259;
    static const double KD = 2.1998873957149437;
    static const double KE = 2.5687780604965515;
    static const double KF = 3.6189656655262064;
    static const double KG = 5.3616789503774974;
    static const double KH = 8.7012952137792582;

    static const double aA = 10.916991008220921;
    static const double aB = 10.226477662739623;
    static const double aC = 7.989784392423307;
    static const double aD = 5.4216606462081334;
    static const double aE = 1.904421638973281;
    static const double aF = 1.1476357111552424;

    static const double bA = 1.6623966587221613;
    static const double bB = 1.8517823702277259;
    static const double bC = 2.0747275727604848;
    static const double bD = 2.3843327281057476;
    static const double bE = 3.0938718630113789;
    static const double bF = 4.4903223079518519;

    if (K < KA)
    {
        if (K >= 0.0)
        {
            scipy_zeros_info solver_stats;
            double mc = brentq(icelk_func, 0.0, 1.0, 2.0e-12, 8.881784197001252e-16, 100, &K, &solver_stats);

            switch (solver_stats.error_num)
            {
            case CONVERGED:
                return 1.0 - 1.0 / mc;
            case SIGNERR:
                feraiseexcept(FE_INVALID);
                return nan("");
            case CONVERR:
                feraiseexcept(FE_INVALID);
                return 1.0 - 1.0 / mc;
            default:
                return nan("");
            }
        }
        else
        {
            feraiseexcept(FE_INVALID);
            return nan("");
        }
    }
    else if (K < KB)
    {
        static const double c[] = {
            19401242299345095.0e-17,
            18194753782220755.0e-17,
            -1148955521044524.0e-17,
            55159476593419.0e-17,
            -2283669657743.0e-17,
            86634143620.0e-17,
            -3105064734.0e-17,
            106925302.0e-17,
            -3574353.0e-17,
            116785.0e-17,
            -3747.0e-17,
            118.0e-17,
        };
        enum
        {
            iend = _Countof(c) - 1
        };

        double x = aA * (K - bA);

        double t[iend + 1];
        t[0] = 1.0;
        t[1] = x;

        for (int i = 2; i <= iend; ++i)
        {
            t[i] = 2.0 * x * t[i - 1] - t[i - 2];
        }

        double sum = 0.0;

        for (int i = iend; i >= 0; --i)
        {
            sum += c[i] * t[i];
        }

        return sum;
    }
    else if (K < KC)
    {
        static const double c[] = {
            48990719448270469.0e-17,
            11715040871065477.0e-17,
            -739311643473559.0e-17,
            34907474599398.0e-17,
            -1406019129690.0e-17,
            51595118298.0e-17,
            -1783762901.0e-17,
            59157015.0e-17,
            -1902485.0e-17,
            59756.0e-17,
            -1842.0e-17,
            56.0e-17,
        };
        enum
        {
            iend = _Countof(c) - 1
        };

        double x = aB * (K - bB);

        double t[iend + 1];
        t[0] = 1.0;
        t[1] = x;

        for (int i = 2; i <= iend; ++i)
        {
            t[i] = 2.0 * x * t[i - 1] - t[i - 2];
        }

        double sum = 0.0;

        for (int i = iend; i >= 0; --i)
        {
            sum += c[i] * t[i];
        }

        return sum;
    }
    else if (K < KD)
    {
        static const double c[] = {
            69305133428229853.0e-17,
            8613490854834244.0e-17,
            -653420708906968.0e-17,
            36436506978491.0e-17,
            -1709827595538.0e-17,
            72527965477.0e-17,
            -2887658427.0e-17,
            110075215.0e-17,
            -4063750.0e-17,
            146392.0e-17,
            -5173.0e-17,
            180.0e-17,
        };
        enum
        {
            iend = _Countof(c) - 1
        };

        double x = aC * (K - bC);

        double t[iend + 1];
        t[0] = 1.0;
        t[1] = x;

        for (int i = 2; i <= iend; ++i)
        {
            t[i] = 2.0 * x * t[i - 1] - t[i - 2];
        }

        double sum = 0.0;

        for (int i = iend; i >= 0; --i)
        {
            sum += c[i] * t[i];
        }

        return sum;
    }
    else if (K < KE)
    {
        static const double c[] = {
            84201552855821970.0e-17,
            6201036207545733.0e-17,
            -648514448920333.0e-17,
            48795524452071.0e-17,
            -3029715923753.0e-17,
            167839458005.0e-17,
            -8670547452.0e-17,
            427595012.0e-17,
            -20387219.0e-17,
            947243.0e-17,
            -43128.0e-17,
            1932.0e-17,
            -85.0e-17,
        };
        enum
        {
            iend = _Countof(c) - 1
        };

        double x = aD * (K - bD);

        double t[iend + 1];
        t[0] = 1.0;
        t[1] = x;

        for (int i = 2; i <= iend; ++i)
        {
            t[i] = 2.0 * x * t[i - 1] - t[i - 2];
        }

        double sum = 0.0;

        for (int i = iend; i >= 0; --i)
        {
            sum += c[i] * t[i];
        }

        return sum;
    }
    else if (K < KF)
    {
        static const double c[] = {
            3372842411096568.0e-15,
            1082799019684612.0e-15,
            -6377916757929.0e-15,
            804553258630.0e-15,
            -72934012596.0e-15,
            5172172149.0e-15,
            -325770457.0e-15,
            21898786.0e-15,
            -1696166.0e-15,
            136104.0e-15,
            -10427.0e-15,
            780.0e-15,
            -60.0e-15,
            5.0e-15,
        };
        enum
        {
            iend = _Countof(c) - 1
        };

        double x = aE * (K - bE);

        double t[iend + 1];
        t[0] = 1.0;
        t[1] = x;

        for (int i = 2; i <= iend; ++i)
        {
            t[i] = 2.0 * x * t[i - 1] - t[i - 2];
        }

        double sum = 0.0;

        for (int i = iend; i >= 0; --i)
        {
            sum += c[i] * t[i];
        }

        return 1.0 - exp(-sum);
    }
    else if (K < KG)
    {
        static const double c[] = {
            620234816387263.0e-14,
            174947955583875.0e-14,
            -226305230130.0e-14,
            50978226844.0e-14,
            -8407174303.0e-14,
            1058249590.0e-14,
            -103492197.0e-14,
            7909234.0e-14,
            -487701.0e-14,
            30082.0e-14,
            -2890.0e-14,
            372.0e-14,
            -42.0e-14,
            4.0e-14,
        };
        enum
        {
            iend = _Countof(c) - 1
        };

        double x = aF * (K - bF);

        double t[iend + 1];
        t[0] = 1.0;
        t[1] = x;

        for (int i = 2; i <= iend; ++i)
        {
            t[i] = 2.0 * x * t[i - 1] - t[i - 2];
        }

        double sum = 0.0;

        for (int i = iend; i >= 0; --i)
        {
            sum += c[i] * t[i];
        }

        return 1.0 - exp(-sum);
    }
    else if (K < KH)
    {
        double kx = K - 6.6363331625866755;
        /* clang-format off */
        double sum = 10.5 + 
            kx * (2.0001414497655 + 
            kx * (-1.277026552e-4 + 
            kx * (7.59811036e-5 + 
            kx * (-3.34231082e-5 + 
            kx * (1.15490151e-5 + 
            kx * (-3.2468658e-6 + 
            kx * (7.573936e-7 + 
            kx * (-1.476109e-7 + 
            kx * (2.383963e-8 + 
            kx * (-3.07883e-9 + 
            kx * (2.83187e-10 + 
            kx * (-9.125e-12 + 
            kx * (-2.572e-12 + 
            kx * (5.45e-13 + 
            kx * (-3.0e-14 + 
            kx * (-2.0e-14 + 
            kx * 6.2e-15))))))))))))))));
        /* clang-format on */
        return 1.0 - exp(-sum);
    }
    else if (K < 19.408121055678468)
    {
        double kx = K - 10.1862944132991;
        /* clang-format off */
        double sum = 17.6 + 
            kx * (2.0000001973566 + 
            kx * (-1.859964e-7 + 
            kx * (1.164242e-7 + 
            kx * (-5.44254e-8 + 
            kx * (2.02555e-8 + 
            kx * (-6.247e-9 + 
            kx * (1.6406e-9 + 
            kx * (-3.74e-10 + 
            kx * (7.5e-11 + 
            kx * -1.3e-11)))))))));
        /* clang-format on */
        return 1.0 - exp(-sum);
    }
    else if (isnan(K))
    {
        return K;
    }
    else
    {
        feraiseexcept(FE_INEXACT);
        return 1.0;
    }
}

static double icele_func(double mc, void *E)
{
    double e = *(double *)E;
    return ceie(mc) - sqrt(mc) * e;
}

double icele(double EP)
{
    static const double EA = 1.5707963267948966;
    static const double EB = 1.4111237670965148;
    static const double EC = 1.2907868683852396;
    static const double ED = 1.1878333138534634;
    static const double EE = 1.1199612906642078;
    static const double EF = 1.0556008611359877;
    static const double EG = 1.0048203911429256;
    static const double EH = 1.000044347637373;

    static const double aA = -12.5256337330469197;
    static const double aB = -16.6200061778109266;
    static const double aC = -19.4262355398589756;
    static const double aD = -29.4672223697113173;
    static const double aE = 2.62838060962536483;
    static const double aF = 0.817880966989805936;

    static const double bA = 1.4909600469457057;
    static const double bB = 1.3509553177408772;
    static const double bC = 1.2393100911193515;
    static const double bD = 1.1538973022588356;
    static const double bE = 2.5090941878294072;
    static const double bF = 4.11222839700788;

    double E = 1.0 + EP;

    if (E > EA)
    {
        if (isinf(E))
        {
            return -E;
        }
        else
        {
            scipy_zeros_info solver_stats;
            double mc = brentq(icele_func, 0.0, 1.0, 2.0e-12, 8.881784197001252e-16, 100, &E, &solver_stats);

            switch (solver_stats.error_num)
            {
            case CONVERGED:
                return 1.0 - 1.0 / mc;
            case SIGNERR:
                feraiseexcept(FE_INVALID);
                return nan("");
            case CONVERR:
                feraiseexcept(FE_INVALID);
                return 1.0 - 1.0 / mc;
            default:
                return nan("");
            }
        }
    }
    else if (E > EB)
    {
        static const double c[] = {
            0.191491596762105577,
            0.187521750862407743,
            -0.003990799269547337,
            -0.000021715474908196,
            -0.000000795715611047,
            -0.000000035290453296,
            -0.000000001771319411,
            -0.000000000096704988,
            -0.000000000005606172,
            -0.000000000000339873,
            -0.00000000000002133,
            -0.000000000000001376,
        };
        enum
        {
            iend = _Countof(c) - 1,
        };

        double x = aA * (E - bA);

        double t[iend + 1];
        t[0] = 1.0;
        t[1] = x;

        for (int i = 2; i <= iend; ++i)
        {
            t[i] = 2.0 * x * t[i - 1] - t[i - 2];
        }

        double sum = 0.0;

        for (int i = iend; i >= 0; --i)
        {
            sum += c[i] * t[i];
        }

        return sum;
    }
    else if (E > EC)
    {
        static const double c[] = {
            0.4969412566602583,
            0.11951688581417314,
            -0.00244061880474704,
            -0.00001685618407984,
            -0.00000063628928178,
            -0.00000002953987502,
            -0.00000000156070789,
            -0.00000000008986478,
            -0.00000000000549838,
            -0.00000000000035191,
            -0.00000000000002332,
            -0.00000000000000159,
        };
        enum
        {
            iend = _Countof(c) - 1,
        };

        double x = aB * (E - bB);

        double t[iend + 1];
        t[0] = 1.0;
        t[1] = x;

        for (int i = 2; i <= iend; ++i)
        {
            t[i] = 2.0 * x * t[i - 1] - t[i - 2];
        }

        double sum = 0.0;

        for (int i = iend; i >= 0; --i)
        {
            sum += c[i] * t[i];
        }

        return sum;
    }
    else if (E > ED)
    {
        static const double c[] = {
            0.70197876373537509,
            0.08602065171825419,
            -0.00197779951710636,
            -0.00002059553190095,
            -0.00000096047713357,
            -0.00000005591458166,
            -0.00000000372015072,
            -0.0000000002700766,
            -0.00000000002084262,
            -0.00000000000168273,
            -0.00000000000014065,
            -0.00000000000001208,
            -0.00000000000000106,
        };
        enum
        {
            iend = _Countof(c) - 1,
        };

        double x = aC * (E - bC);

        double t[iend + 1];
        t[0] = 1.0;
        t[1] = x;

        for (int i = 2; i <= iend; ++i)
        {
            t[i] = 2.0 * x * t[i - 1] - t[i - 2];
        }

        double sum = 0.0;

        for (int i = iend; i >= 0; --i)
        {
            sum += c[i] * t[i];
        }

        return sum;
    }
    else if (E > EE)
    {
        static const double c[] = {
            0.83448854563061028,
            0.04751257054900759,
            -0.00098795130570153,
            -0.00001253506060813,
            -0.00000059189763242,
            -0.00000003530719901,
            -0.00000000241289912,
            -0.00000000018000725,
            -0.00000000001427478,
            -0.00000000000118417,
            -0.00000000000010167,
            -0.00000000000000897,
            -0.00000000000000081,
        };
        enum
        {
            iend = _Countof(c) - 1,
        };

        double x = aD * (E - bD);

        double t[iend + 1];
        t[0] = 1.0;
        t[1] = x;

        for (int i = 2; i <= iend; ++i)
        {
            t[i] = 2.0 * x * t[i - 1] - t[i - 2];
        }

        double sum = 0.0;

        for (int i = iend; i >= 0; --i)
        {
            sum += c[i] * t[i];
        }

        return sum;
    }
    else if (E > EF)
    {
        static const double c[] = {
            2.621927976267172,
            0.480434270533396,
            -0.002421516965161,
            0.000060368511503,
            -0.000001103238537,
            0.000000004267375,
            0.00000000065424,
            -0.000000000030728,
            0.000000000000745,
            -0.000000000000005,
        };
        enum
        {
            iend = _Countof(c) - 1,
        };

        double r = -log(EP);
        double x = aE * (r - bE);

        double t[iend + 1];
        t[0] = 1.0;
        t[1] = x;

        for (int i = 2; i <= iend; ++i)
        {
            t[i] = 2.0 * x * t[i - 1] - t[i - 2];
        }

        double sum = 0.0;

        for (int i = iend; i >= 0; --i)
        {
            sum += c[i] * t[i];
        }
        return 1.0 - exp(-sum);
    }
    else if (E > EG)
    {
        static const double c[] = {
            4.563857924247435,
            1.448978332531438,
            -0.013787853000146,
            0.001017879586053,
            -0.000069972247182,
            0.000003797035826,
            -0.000000100760933,
            -0.000000008981925,
            0.000000001749497,
            -0.000000000170996,
            0.000000000011351,
            -0.0000000000004,
            -0.000000000000022,
        };
        enum
        {
            iend = _Countof(c) - 1,
        };

        double r = -log(EP);
        double x = aF * (r - bF);

        double t[iend + 1];
        t[0] = 1.0;
        t[1] = x;

        for (int i = 2; i <= iend; ++i)
        {
            t[i] = 2.0 * x * t[i - 1] - t[i - 2];
        }

        double sum = 0.0;

        for (int i = iend; i >= 0; --i)
        {
            sum += c[i] * t[i];
        }

        return 1.0 - exp(-sum);
    }
    else if (E > EH)
    {
        double ex = -log(EP) - 7.106602166777170551;
        /* clang-format off */
        double sum = 8.0 + 
            ex * (1.1138558445808138764 + 
            ex * (-7.160737009252904412e-3 + 
            ex * (6.14816637761486272e-4 + 
            ex * (-5.7949395860221721e-5 + 
            ex * (5.45152607032739e-6 + 
            ex * (-4.80075731204772e-7 + 
            ex * (3.7093529426652e-8 + 
            ex * (-2.222894130707e-9 + 
            ex * (5.37796726e-11 + 
            ex * (1.071061734e-11 + 
            ex * (-2.33974846e-12 + 
            ex * (3.0858084e-13 + 
            ex * (-3.239999e-14 + 
            ex * (2.84089e-15 + 
            ex * (-2.0071e-16 + 
            ex * (9.0e-18 + 
            ex * (3.0e-19 + 
            ex * -2.0e-19)))))))))))))))));
        /* clang-format on */
        return 1.0 - exp(-sum);
    }
    else
    {
        double ex = -log(EP) - 13.566548324200306727;
        /* clang-format off */
        double sum = 15.0 + 
            ex * (1.0634010131365226997 + 
            ex * (-2.137214443579717823e-3 + 
            ex * (9.8905527968533563e-5 + 
            ex * (-5.240426861207633e-6 + 
            ex * (2.99271455626484e-7 + 
            ex * (-1.7869634729734e-8 + 
            ex * (1.092804281815e-9 + 
            ex * (-6.7254295377e-11 + 
            ex * (4.096019951e-12 + 
            ex * (-2.42723757e-13 + 
            ex * (1.3740808e-14 + 
            ex * (-7.2627e-16 + 
            ex * (3.455e-17 + 
            ex * (-1.36e-18 + 
            ex * (3.0e-20 + 
            ex * (-1.4e-21 + 
            ex * 7.e-23))))))))))))))));
        /* clang-format on */
        return 1.0 - exp(-sum);
    }
}

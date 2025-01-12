#ifndef YNOGK_CXX_PARTICLE_HXX
#define YNOGK_CXX_PARTICLE_HXX

#include <vector>

#include "particle.h"

namespace ynogk
{
    class Particle
    {
    public:
        Particle(double a_spin, double robs, double muobs, double sinobs, double scal, const double *velocity);

        Particle(const Particle &);
        Particle(Particle &&);
        ~Particle();

        void lambdaq(double alpha, double beta);
        void lambdaq(double pr, double ptheta, double pphi);
        void metric();
        void center_of_image();
        int mutp();
        int rtp();
        double mu2p(double mu, int t1, int t2);
        double r2p(double r, int t1, int t2);
        double mucos(double pem);
        double radius(double pem);
        double phi(double pem);
        double pemdisk(double mu, double rout, double rin);
        double pemdisk_all(double mu, double rout, double rin);
        double p_total();
        void ynogk(double pem);

        auto operator->()
        {
            return &this->pt;
        }

    private:
        ptcl pt;
        class Impl;
        std::vector<Impl> pimpl;
    };
}

#endif

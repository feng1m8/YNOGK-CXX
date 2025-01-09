#include "impl.hxx"

namespace ynogk
{
    Particle::Particle(double a_spin, double robs, double muobs, double sinobs, double scal, double *velocity) : pimpl(1, a_spin)
    {
        particle_construct(&this->pt, a_spin, robs, muobs, sinobs, scal, velocity);
    }

    Particle::Particle(const Particle &) = default;
    Particle::Particle(Particle &&) = default;
    Particle::~Particle() = default;

    void Particle::lambdaq(double alpha, double beta)
    {
        this->pt.alpha = alpha;
        this->pt.beta = beta;
        ::lambdaq(&this->pt);
    }

    void Particle::lambdaq(double pr, double ptheta, double pphi)
    {
        ini_direction2lamdaq(&this->pt, pr, ptheta, pphi);
    }

    void Particle::metric()
    {
        metricg(&this->pt);
    }

    void Particle::center_of_image()
    {
        ::center_of_image(&this->pt);
    }

    int Particle::mutp()
    {
        return ::mutp(&this->pt);
    }

    int Particle::rtp()
    {
        return radiustp(&this->pt);
    }

    double Particle::mu2p(double mu, int t1, int t2)
    {
        return ::mu2p(&this->pt, mu, t1, t2);
    }

    double Particle::r2p(double r, int t1, int t2)
    {
        return this->pimpl.front().radius.coordinates.r2p(&this->pt, r, t1, t2);
    }

    double Particle::mucos(double pem)
    {
        return this->pimpl.front().radius.coordinates.mucos(&this->pt, pem);
    }

    double Particle::radius(double pem)
    {
        return this->pimpl.front().radius.coordinates.radius(&this->pt, pem);
    }

    double Particle::phi(double pem)
    {
        return this->pimpl.front().phi(&this->pt, pem);
    }

    double Particle::pemdisk(double mu, double rout, double rin)
    {
        return this->pimpl.front().radius.coordinates.Pemdisk(&this->pt, mu, rout, rin, &this->pt.radius);
    }

    double Particle::pemdisk_all(double mu, double rout, double rin)
    {
        return this->pimpl.front().radius.coordinates.Pemdisk_all(&this->pt, mu, rout, rin, &this->pt.radius);
    }

    double Particle::p_total()
    {
        return this->pimpl.front().radius.coordinates.p_total(&this->pt);
    }

    void Particle::ynogk(double pem)
    {
        this->pimpl.front().YNOGKC(&this->pt, pem);
    }
}

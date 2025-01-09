cimport cython
cimport cython.operator as operator
from libcpp cimport cmath


cdef extern from 'ynogk_cxx/particle.hxx':
    'static auto &operator*(ynogk::Particle &self) { return *self.operator->(); }'

    double c_rms 'rms' (double) nogil
    int metricgij(double, double, double, double, double *, double *, double *, double *, double *) nogil

    ctypedef struct metrics:
        double somiga
        double expnu
        double exppsi
        double expmu1
        double expmu2

    ctypedef struct c_ptcl 'ptcl':
        double robs
        double muobs
        double sinobs
        double[4] velocity_ini
        double scal
        double a_spin
        double rhorizon
        double[5] f1234
        double alpha
        double beta
        double lamda 'lambda'
        double q
        double alphac
        double betac
        double mu_tp1
        double mu_tp2
        unsigned int mu_reals
        int mobseqmtp
        double r_tp1
        double r_tp2
        unsigned int r_reals
        int robs_eq_rtp
        int indrhorizon
        int cases
        double r_p
        double mu_p
        double sin_p
        double phi_p
        double time_p
        double sigma_p
        double sign_pr_p
        double sign_pth_p
        double radius
        double sign_pr
        double sign_pth
        metrics mt

    cdef cppclass Particle 'ynogk::Particle':
        Particle(double, double, double, double, double, double *) except +
        void lambdaq(double, double)
        void lambdaq(double, double, double)
        void metric()
        void center_of_image()
        int mutp()
        int rtp()
        double mu2p(double, int, int)
        double r2p(double, int, int)
        double mucos(double)
        double radius(double)
        double phi(double)
        double pemdisk(double, double, double)
        double pemdisk_all(double, double, double)
        double p_total()
        void ynogk(double)
        c_ptcl operator*()


@cython.ufunc
cdef double rms(double a_spin) nogil:
    return c_rms(a_spin)


@cython.ufunc
cdef double rhorizon(double a_spin) nogil:
    return 1.0 + cmath.sqrt(1.0 - a_spin)


@cython.ufunc
cdef (double, double, double, double, double) metric(double robs, double muobs, double sinobs, double a_spin) noexcept nogil:
    cdef double omega, expnu, exppsi, expmu1, expmu2
    metricgij(robs, muobs, sinobs, a_spin, &omega, &expnu, &exppsi, &expmu1, &expmu2)
    return omega, expnu, exppsi, expmu1, expmu2


@cython.cpp_locals(True)
cdef class ptcl:
    cdef Particle pt
    
    def __init__(self, double a_spin, double radius, double muobs, double sinobs, velocity=[0.0, 0.0, 0.0], double scal=1.0):
        cdef double[3] velo = [velocity[0], velocity[1], velocity[2]]
        self.pt = Particle(a_spin, radius, muobs, sinobs, scal, velo)

    def lambdaq(self, double pr, double ptheta, pphi=None):
        if pphi is None:
            self.pt.lambdaq(pr, ptheta)
        else:
            self.pt.lambdaq(pr, ptheta, pphi)

    def metric(self):
        self.pt.metric()

    def center_of_image(self):
        self.pt.center_of_image()

    def mutp(self):
        return self.pt.mutp()    

    def rtp(self):
        return self.pt.rtp()    

    def mu2p(self, double mu, int t1, int t2):
        return self.pt.mu2p(mu, t1, t2)    

    def r2p(self, double r, int t1, int t2):
        return self.pt.r2p(r, t1, t2)    

    def mucos(self, double pem):
        operator.dereference(self.pt).mu_p = self.pt.mucos(pem)
        operator.dereference(self.pt).sin_p = cmath.sqrt(1.0 - operator.dereference(self.pt).mu_p * operator.dereference(self.pt).mu_p)
        return operator.dereference(self.pt).mu_p

    def radius(self, double pem):
        return self.pt.radius(pem)

    def phi(self, double pem):
        self.pt.phi(pem)
        operator.dereference(self.pt).radius = operator.dereference(self.pt).r_p
        operator.dereference(self.pt).sign_pr = operator.dereference(self.pt).sign_pr_p
        operator.dereference(self.pt).sign_pth = operator.dereference(self.pt).sign_pth_p
        return operator.dereference(self.pt).phi_p

    def pemdisk(self, double mu, double rout, double rin):
        return self.pt.pemdisk(mu, rout, rin)

    def pemdisk_all(self, double mu, double rout, double rin):
        return self.pt.pemdisk_all(mu, rout, rin)
        
    def p_total(self):
        return self.pt.p_total()

    def ynogk(self, double pem):
        self.pt.ynogk(pem)
        operator.dereference(self.pt).radius = operator.dereference(self.pt).r_p
        operator.dereference(self.pt).sign_pr = operator.dereference(self.pt).sign_pr_p
        operator.dereference(self.pt).sign_pth = operator.dereference(self.pt).sign_pth_p

    @property
    def robs(self):
        return operator.dereference(self.pt).robs

    @property
    def muobs(self):
        return operator.dereference(self.pt).muobs

    @property
    def sinobs(self):
        return operator.dereference(self.pt).sinobs

    @property
    def velocity(self):
        return [operator.dereference(self.pt).velocity_ini[1], operator.dereference(self.pt).velocity_ini[2], operator.dereference(self.pt).velocity_ini[3]]

    @property
    def scal(self):
        return operator.dereference(self.pt).scal

    @property
    def a_spin(self):
        return operator.dereference(self.pt).a_spin

    @property
    def rhorizon(self):
        return operator.dereference(self.pt).rhorizon

    @property
    def f1234(self):
        return [operator.dereference(self.pt).f1234[4], operator.dereference(self.pt).f1234[1], operator.dereference(self.pt).f1234[2], operator.dereference(self.pt).f1234[3]]

    @property
    def alpha(self):
        return operator.dereference(self.pt).alpha

    @property
    def beta(self):
        return operator.dereference(self.pt).beta

    @property
    def lamda(self):
        return operator.dereference(self.pt).lamda

    @property
    def q(self):
        return operator.dereference(self.pt).q

    @property
    def alphac(self):
        return operator.dereference(self.pt).alphac

    @property
    def betac(self):
        return operator.dereference(self.pt).betac

    @property
    def mu_tp1(self):
        return operator.dereference(self.pt).mu_tp1

    @property
    def mu_tp2(self):
        return operator.dereference(self.pt).mu_tp2

    @property
    def mu_reals(self):
        return operator.dereference(self.pt).mu_reals

    @property
    def mu_is_mutp(self):
        if operator.dereference(self.pt).mobseqmtp == 1:
            return True
        else:
            return False

    @property
    def r_tp1(self):
        return operator.dereference(self.pt).r_tp1

    @property
    def r_tp2(self):
        return operator.dereference(self.pt).r_tp2

    @property
    def r_reals(self):
        return operator.dereference(self.pt).r_reals

    @property
    def r_is_rtp(self):
        if operator.dereference(self.pt).robs_eq_rtp == 1:
            return True
        else:
            return False

    @property
    def r_to_horizon(self):
        if operator.dereference(self.pt).indrhorizon == 1:
            return True
        else:
            return False

    @property
    def r_to_inf(self):
        if operator.dereference(self.pt).cases == 1:
            return True
        else:
            return False

    @property
    def r_p(self):
        return operator.dereference(self.pt).radius

    @property
    def mu_p(self):
        return operator.dereference(self.pt).mu_p

    @property
    def sin_p(self):
        return operator.dereference(self.pt).sin_p

    @property
    def phi_p(self):
        return operator.dereference(self.pt).phi_p

    @property
    def time_p(self):
        return operator.dereference(self.pt).time_p

    @property
    def sigma_p(self):
        return operator.dereference(self.pt).sigma_p

    @property
    def sign_pr_p(self):
        return operator.dereference(self.pt).sign_pr

    @property
    def sign_ptheta_p(self):
        return operator.dereference(self.pt).sign_pth

    @property
    def omega(self):
        return operator.dereference(self.pt).mt.somiga

    @property
    def expnu(self):
        return operator.dereference(self.pt).mt.expnu

    @property
    def exppsi(self):
        return operator.dereference(self.pt).mt.exppsi

    @property
    def expmu1(self):
        return operator.dereference(self.pt).mt.expmu1

    @property
    def expmu2(self):
        return operator.dereference(self.pt).mt.expmu2

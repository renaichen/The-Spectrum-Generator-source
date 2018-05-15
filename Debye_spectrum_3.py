import numpy as np
import random


class Generator(object):
    """
    Class which implements the generation of the random part
    and damping part (with memory kernel) for displacement of
    anchor bath particle y_l according to the paper: J. Chem. Phys. 69, 2525 (1978)
    """
    def __init__(self, n, omegaD,
                 mass=1, temperature=1, dt=0.001, t_num=1000, Ntraj=1, sampling_rate=2):

        # set the parameters and constants
        self.n = n
        self.omegaD = omegaD
        self.kB = 0.001987   #kcal/(mol K)
        # self.kB = 0.008314  # kJ/(mol K)
        self.mass = mass
        self.temperature = temperature
        self.dt = dt
        self.Ntraj = Ntraj

        # Initialize the coefficients
        self._ap = np.zeros(self.n)
        self._bp = np.zeros(self.n)
        self._sigma = 1/(2*np.sin((3*np.pi)/(2.0*self.n)))
        self._coe_rho = (self.kB * self.temperature / self.mass) * (
                            2 ** (2 * self.n - 1) * np.sin(
                                 3 * np.pi / (2 * self.n))) / self.omegaD ** 3  # for the coefficient of random number
        self.sspold = np.zeros(self.n)
        self.scpold = np.zeros(self.n)
        self.ssp = np.zeros(self.n)
        self.scp = np.zeros(self.n)

        # set dynamical variables
        self.t_num = t_num
        self.sampling = t_num // sampling_rate
        self.t = np.linspace(0, self.dt * self.t_num, self.t_num)
        self.R_length = np.zeros(self.t_num)
        self.R_sampling = np.zeros(self.sampling)


    def compute_ap(self):
        for p in range(self.n):
            self._ap[p] = np.cos((2 * p + 1) * np.pi / (2.0 * self.n))

    def compute_bp(self):
        for p in range(self.n):
            self._bp[p] = np.sin((2 * p + 1) * np.pi / (2.0 * self.n))

    def compute_coe_rho(self):
        for i in range(1, self.n):
            self._coe_rho *= np.sin(i * np.pi / (2 * self.n)) ** 2

    def get_coefficients(self):
        self.compute_ap()
        self.compute_bp()
        self.compute_coe_rho()
        # return self._coe_rho

    # def damp_evolve_euler(self, bp, ap, ft):
    #     """"
    #     this function is to calculate the ssp and scp in one-time 'damping' part,
    #     that is u(t) in the paper, for n=even case
    #     """
    #     ssp = 1.0
    #     ssptemp = 0.0
    #     scp = 1.0
    #     scptemp = 0.0
    #     tolerance = 1e-10
    #     while abs((ssp - ssptemp) / ssp) > tolerance or abs((scp - scptemp) / scp) > tolerance:
    #         ssptemp = ssp
    #         scptemp = scp
    #         ssp = ssptemp + (-bp * self.omegaD * ssptemp + ap * self.omegaD * scptemp) * self.dt
    #         scp = scptemp + (-bp * self.omegaD * scptemp - ap * self.omegaD * ssptemp + ft) * self.dt
    #         # print abs((ssp - ssptemp) / ssp)
    #     return ssp, scp,
    #
    # def damp_getter(self, ft):   # Only even case is taken into account for now
    #     self.get_coefficients()
    #     ssp = np.zeros(self.n//2)
    #     scp = np.zeros(self.n//2)
    #     ut = 0.0
    #     for i in range(self.n//2):
    #         ssp[i], scp[i] = self.damp_evolve_euler(self._bp[i], self._ap[i], ft)
    #         ut += 2 * self._ap[i] * self._bp[i] * ssp[i] + (self._bp[i]**2 - self._ap[i]**2) * scp[i]
    #     return ut * (1.0/(self.mass * self.omegaD * self._sigma))

    def damp_getter(self, ft):   # Only even case is taken into account for now
        self.compute_ap()
        self.compute_bp()

        # ut = 0.0
        self.ssp = self.sspold + (-self._bp * self.omegaD * self.sspold +
                                    self._ap * self.omegaD * self.scpold) * self.dt
        self.scp = self.scpold + (-self._bp * self.omegaD * self.scpold -
                                    self._ap * self.omegaD * self.sspold + ft) * self.dt
        # ut += 2 * self._ap[i] * self._bp[i] * ssp[i] + (self._bp[i]**2 - self._ap[i]**2) * scp[i]
        uttwice = 2 * self._ap * self._bp * self.ssp + (self._bp**2 - self._ap**2) * self.scp
        ut = np.sum(uttwice[: self.n//2])

        self.sspold = self.ssp
        self.scpold = self.scp
        # print self.sspold
        return ut * (1.0/(self.mass * self.omegaD * self._sigma))

    def random_evolve_euler(self, x, v):
        tstep = 0
        while tstep < (self.t_num-1):
            # x[tstep, self.n / 2] = self._coe_rho * random.gauss(0, 1)
            x[tstep, self.n / 2] = np.sqrt(self._coe_rho) * random.gauss(0, 1)

            # This x[..., n/2] elements, they are the original Gaussian random numberss,
            # so the following x[..., i] should be treated differently than x[..., i-1]

            for i in range(self.n / 2, 0, -1):
                # x[tstep + 1, i - 1] = x[tstep, i - 1] + v[tstep, i - 1] * self.dt
                x[tstep + 1, i - 1] = x[tstep, i - 1] + v[tstep, i - 1] * self.dt
                # x[tstep + 1, i - 1] = np.exp(-1. / 5000 * self.t[tstep + 1]) * x[tstep + 1, i - 1]
                # v[tstep + 1, i - 1] = v[tstep, i - 1] + self.omegaD ** 2 * \
                #                         (x[tstep, i] - x[tstep, i - 1] -
                #                             2 * self._bp[i - 1] / self.omegaD * v[tstep, i - 1]) * self.dt

                #---------------------------------------------------
                # v[tstep + 1, i - 1] = v[tstep, i - 1] + self.omegaD ** 2 * \
                #                                         (x[tstep, i] / np.sqrt(self.dt) - x[tstep, i - 1] -
                #                                          2 * self._bp[i - 1] / self.omegaD * v[tstep, i - 1]) * self.dt
                ##---------------------------------------------------

                ##----- One must use the following condition (not the one above)
                # to do the none-gaussian velocities deterministically, otherwise it blows up
                if i == self.n / 2:
                    v[tstep + 1, i - 1] = v[tstep, i - 1] + self.omegaD ** 2 * \
                                                        (x[tstep, i]/np.sqrt(self.dt) - x[tstep, i - 1] -
                                                         2 * self._bp[i - 1] / self.omegaD * v[tstep, i - 1]) * self.dt
                else:
                    v[tstep + 1, i - 1] = v[tstep, i - 1] + self.omegaD ** 2 * \
                                                        (x[tstep, i] - x[tstep, i - 1] -
                                                        2 * self._bp[i - 1] / self.omegaD * v[tstep, i - 1]) * self.dt


                # v[tstep + 1, i - 1] = np.exp(-1./6000*self.t[tstep+1])*v[tstep+1, i - 1]

            tstep += 1

    def random_evolve_vv(self, x, v):
        tstep = 0
        a = np.zeros((self.t_num, self.n / 2))
        while tstep < (self.t_num-1):
            # x[tstep, self.n / 2] = self._coe_rho * random.gauss(0, 1)
            x[tstep, self.n / 2] = np.sqrt(self._coe_rho) * random.gauss(0, 1)

            # This x[..., n/2] elements, they are the original Gaussian random numberss,
            # so the following x[..., i] should be treated differently than x[..., i-1]

            for i in range (self.n/2, 0, -1):
                a[tstep, i - 1] = - self.omegaD ** 2 * x[tstep, i - 1]
                x[tstep + 1, i - 1] = x[tstep, i - 1] + v[tstep, i - 1] * self.dt + 0.5 * a[tstep, i - 1] * self.dt ** 2
                a[tstep + 1, i - 1] = - self.omegaD ** 2 * x[tstep + 1, i - 1]

                ##------------------------------------------------------------
                # v[tstep + 1, i - 1] = v[tstep, i - 1] + 0.5 * (a[tstep, i - 1] + a[tstep + 1, i - 1]) * self.dt \
                #                       - 2 * self._bp[i - 1] * self.omegaD * v[tstep, i - 1] * self.dt \
                #                       + self.omegaD ** 2 * x[tstep, i] * np.sqrt(self.dt)
                ##---------------------------------------------------

                ##----- One must use the following condition (not the one above)
                # to do the none-gaussian velocities deterministically, otherwise it blows up
                if i == self.n / 2:
                    v[tstep + 1, i - 1] = v[tstep, i - 1] + 0.5 * (a[tstep, i - 1] + a[tstep + 1, i - 1]) * self.dt \
                                          - 2 * self._bp[i - 1] * self.omegaD * v[tstep, i - 1] * self.dt \
                                          + self.omegaD ** 2 * x[tstep, i] * np.sqrt(self.dt)
                else:
                    v[tstep + 1, i - 1] = v[tstep, i - 1] + 0.5 * (a[tstep, i - 1] + a[tstep + 1, i - 1]) * self.dt \
                                          - 2 * self._bp[i - 1] * self.omegaD * v[tstep, i - 1] * self.dt \
                                          + self.omegaD ** 2 * x[tstep, i] * self.dt

            tstep += 1

    def random_mult_traj(self):
        self.compute_bp()
        self.compute_coe_rho()
        sampling = self.sampling
        self.R_traj = np.zeros((sampling, self.Ntraj))
        if self.n % 2 == 0:     # it is even
            traj = 0
            while traj<self.Ntraj:
                """
                zz is used to represent the Z matrix that is decomposed from the set of 2nd-DEs,
                dzdt is just the derivative elements of that
                """
                zz = np.zeros((self.t_num, self.n / 2 + 1))
                dzdt = np.zeros((self.t_num, self.n / 2 + 1))
                zz[0, :] = np.ones(self.n / 2 + 1)

                self.random_evolve_euler(zz, dzdt)
                # self.random_evolve_vv(zz, dzdt)

                self.R_traj[:, traj] = zz[self.t_num - sampling :, 0]
                self.R_length[:] += zz[:, 0]
                print "trajectory for better random position sampling #: ", traj
                traj += 1
            self.R_sampling[:] = self.R_length[self.t_num - sampling :]/self.Ntraj
            self.Raver = np.sum(self.R_sampling)/(sampling)
        else:    # it is odd
            pass

    def give_me_random_series(self, dt):
        # Need to feed in time step length for the MD simulation, in unit of ps
        division = np.true_divide(dt, self.dt)  # This numpy division operator is better than // or %
        new_points = int(division)
        if dt < self.dt:
            print "Cannot generate random series with " \
                  "step smaller than our iteraction step... " \
                  "please use a smaller time step for the Generator"
            exit()
        elif division != new_points:
            print "the MD time step cannot be divided by Generator " \
                  "time step completely, use a smaller Generator step maybe..."
            exit()
        else:
            self.random_mult_traj()
            rand_out = np.array([])
            i = 0
            while i * new_points < len(self.R_sampling):
                rand_out = np.append(rand_out, self.R_sampling[i * new_points])
                i += 1
            return rand_out


###===================Generator class ends========================


##====================functions defined outside the class==========================
def generate_autocorrelation(seq):
    n = len(seq)
    if n % 2 == 0:
        length = n // 2
    else:
        length = n // 2 + 1
        seq = np.append(seq, 0.)

    correlation = np.zeros(length)
    for i in range(length):
        seq_shift = np.roll(seq, i)
        seq_shift[:i] = np.zeros(i)
        seq_shift[i+length:] = np.zeros(length-i)
        correlation[i] = np.dot(seq, seq_shift) / float(length)

    return correlation

def fourier_transform(seq, deltat):
    N = len(seq)
    # xf = np.linspace(0, 1/(2*deltat), N//2)
    xf = np.linspace(0, np.pi / deltat, N // 2)
    # The normalization for omega is 2pi/deltat,
    # because we use half of the overall points, so we h
    yf = fft(seq)

    # return xf, 2.0 / N * np.abs(yf[0: N // 2])
    # return xf, np.abs(yf[0: N // 2]) / np.abs(yf[0])
    return xf, yf[0: N // 2].real / yf[0].real

#=====================================================================


##================== main ========================================


import sys
import os

import numpy as np
import time
import random
from multiprocessing import Pool

#------------------------------------------------------

def diatomic_traj(n):

    K1 = 0.0

    powerL = 0.0
    powerR = 0.0

    f1new = 0.0
    fL = 0.0
    fR = 0.0
    xL = 46.
    x1new = 50.
    xR = 54.
    v1new = 0.0

    tstep = 0
    while tstep < (tsize-1):
#
        f1old = f1new
        fLold = fL
        fRold = fR
        x1old = x1new
        v1old = v1new

        xiL = random.gauss(0, 1)
        xiR = random.gauss(0, 1)

# #-----------EOM integrator: using the velocity verlet algorithm-----------------------
        x1new = x1old + v1old*dt + (0.5/m1)*f1old*dt**2
        # fL = - AL * alphaL * np.exp(-alphaL * (x1[tstep + 1] - xL[tstep + 1]))
        # fR = AR * alphaR * np.exp(-alphaR * (xR[tstep + 1] - x2[tstep + 1]))

        fL = k1l * (x1new - halfd - xL - x1l)
        fR = - k1r * (xR - x1new - halfd - x1r)

        if x1new < xL or x1new > xR:
            f1 = open('wrong-dia-' + time.strftime('-%m-%d-%H%M%S.txt'), 'w')
            print >> f1, 'error: position disorder, exiting...', omega1, \
                xL, x1new, xR 
            f1.close()
            break

        f1new = -(fL + fR)
#
        v1new = v1old + 0.5*((f1old+f1new)/m1) * dt - \
                gammaL * v1old * dt + np.sqrt(2 * kB * gammaL * TL * 
                            dt / m1) * xiL - \
                gammaR * v1old * dt + np.sqrt(2 * kB * gammaR * TR *
                            dt / m1) * xiR
# #----------------------------------------------------------------------------------
        if tstep > (tsize / 2 - 1):
            powerL += - gammaL * m1 * v1old * 0.5 * (v1old + \
                    v1new) + np.sqrt(2 * kB * gammaL * TL * m1 / dt) * \
                    xiL * 0.5 * (v1old + v1new)
            powerR += - gammaR * m1 * v1old * 0.5 * (v1old + \
                    v1new) + np.sqrt(2 * kB * gammaR * TR * m1 / dt) * \
                    xiR * 0.5 * (v1old + v1new)

            K1 += 0.5 * m1 * v1old ** 2

        tstep += 1

    return K1, powerL, powerR

if __name__ == '__main__':

    SLOTS = int(os.getenv('NSLOTS')) # Get the NSLOTS environment variable provided by the scheduler

    start_time = time.time()

    traj = 3000

    tBegin = 0.
    tEnd = 10000
    dt = 0.01
    tArray = np.arange(tBegin, tEnd, dt)
    tsize = len(tArray)
    halftsize = tsize / 2
    K1traj = np.zeros(traj)
    powerLtraj = np.zeros(traj)
    powerRtraj = np.zeros(traj)

    m1 = float(sys.argv[1])
    mL = 32.
    mR = 32.
    omegaD = 1.
    omegaL = 1.
    omegaR = 1.
    # omega1 = float(sys.argv[1])
    k1l = m1 * omegaL**2
    k1r = m1 * omegaR**2
    # x10 = 50
    x012 = 2.0
    halfd = 1.0
    x1l = 3.0
    x1r = 3.0
    AL = 1 * 1e6  #  1eV = 23kcal/mol
    alphaL = 5e0
    AR = 1 * 1e6  #  1eV = 23kcal/mol
    alphaR = 5e0
    kB = 0.001987
    TL = 270
    TR = 200
    gammaL = 0.1
    gammaR = 0.1
    # gammaL = 3 * omegaL**4 * np.pi / (2 * m1 * mL * omegaD**3)
    # gammaR = 3 * omegaR**4 * np.pi / (2 * m2 * mR * omegaD**3)


    p = Pool(processes=SLOTS)# pass the number of core to the Pool so that I know how many cores I can use.


## Two ways (with index i or using enumerate by remember using tuple() sign)
## seem all good for indexing the trajectory number, but at the safe side
## using i++ is chosen here.
    i = 0
    for K1, powerL, powerR in p.map(diatomic_traj, range(traj)):
    # for i, (K1, K2, K, powerL, power12, powerR) in enumerate(p.map(diatomic_traj, range(traj))):
        K1traj[i] = K1 / halftsize
        powerLtraj[i] = powerL / halftsize
        powerRtraj[i] = powerR / halftsize
        i += 1

    p.close()
    p.join()

    ##----------
    T1aver = np.mean(K1traj) * 2 / kB
    T1std = np.std(K1traj) * 2 / kB
    print 'T1 = ', T1aver, T1std

    JLaver = np.mean(powerLtraj)
    JLstd = np.std(powerLtraj)
    print 'heatL = ', JLaver, JLstd
    JRaver = np.mean(powerRtraj)
    JRstd = np.std(powerRtraj)
    print 'heatR = ', JRaver, JRstd

    run_time = time.time() - start_time
    print 'run time is: ', run_time / 60.

##-----------write-data-out---------
filename = 'center-' + str(traj) + time.strftime('-%m-%d-%H%M%S.txt')
with open(filename, "w") as f:
    f.write("time spent in minutes: %f\n" %(run_time/60))
    # f.write("AL = %f, alphaL = %f\n" %(AL, alphaL))
    f.write("mass = %f\n" %(m1))
    # f.write("equilibrium length (bond length): %f\n" %(x012))
    f.write("trajectory number: %d\n" %(traj))
    f.write("time_step: %f\n" %(dt))
    f.write("number of steps: %d\n" %(tsize/2))
    f.write("TL = %d, TR = %d\n" %(TL, TR))
    f.write("T1 = %f, T1std = %f\n" %(T1aver, T1std))
    f.write("JL = %f, STDJL = %f\n" %(JLaver, JLstd))
    f.write("JR = %f, STDJR = %f\n" %(JRaver, JRstd))

# filename2 = 'heatflux-diatomic-' + str(m1) + time.strftime('-%m-%d-%H%M%S.txt')
# np.savetxt(filename2, np.c_[PsteadyL, Psteady12, PsteadyR])


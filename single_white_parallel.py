import sys
import os

import numpy as np
import time
import random
from multiprocessing import Pool

#------------------------------------------------------

def single_traj(n):

    x1 = np.zeros(tsize)
    xL = np.zeros(tsize)
    xR = np.zeros(tsize)

    v1 = np.zeros(tsize)

    U = np.zeros(tsize)
    K = np.zeros(tsize)
    K1 = np.zeros(tsize)
    UintL = np.zeros(tsize)
    UintR = np.zeros(tsize)

    fLt = np.zeros(tsize)
    fRt = np.zeros(tsize)
    powerL = np.zeros(tsize)
    powerR = np.zeros(tsize)
    powerLsq = np.zeros(tsize)
    powerRsq = np.zeros(tsize)

    xL[0] = 46.
    x1[0] = 50.
    xR[0] = 54.
    halfd = 1.

    v1[0] = 0.0

    f1new = 0.0
    fL = 0.0
    fR = 0.0

    tstep = 0
    while tstep < (tsize-1):
#
        fLold = fL
        fRold = fR
        f1old = f1new

        xiL = random.gauss(0, 1)
        xiR = random.gauss(0, 1)

# #-----------EOM integrator: using the velocity verlet algorithm-----------------------
        x1[tstep+1] = x1[tstep] + v1[tstep]*dt + (0.5/m1)*f1old*dt**2
        xL[tstep+1] = xL[0]
        xR[tstep+1] = xR[0]

        fL = - AL * alphaL * np.exp(-alphaL * (x1[tstep + 1] - halfd - xL[tstep + 1]))
        fR = AR * alphaR * np.exp(-alphaR * (xR[tstep + 1] - halfd - x1[tstep + 1]))

        if x1[tstep+1] < xL[tstep+1] or x1[tstep+1] > xR[tstep+1]:
            f1 = open('wrong-single-' + str(m1) + time.strftime('-%m-%d-%H%M%S.txt'), 'w')
            print >> f1, 'error: position disorder, exiting...', \
                xL[tstep+1], x1[tstep+1], xR[tstep+1] 
            f1.close()
            break

        fLt[tstep] = fL
        fRt[tstep] = fR

        f1new = - (fL + fR)
#
        v1[tstep+1] = v1[tstep] + 0.5*((f1old+f1new)/m1) * dt - \
                        gammaL * v1[tstep] * dt + np.sqrt(2 * kB * gammaL * TL
                                * dt / m1) * xiL - \
                        gammaR * v1[tstep] * dt + np.sqrt(2 * kB * gammaR * TR 
                                * dt / m1) * xiR
# #----------------------------------------------------------------------------------
        UintL[tstep] = AL * np.exp(-alphaL * (x1[tstep] - halfd - xL[tstep]))
        UintR[tstep] = AR * np.exp(-alphaR * (xR[tstep] - halfd - x1[tstep]))
        U[tstep] += UintL[tstep]
        U[tstep] += UintR[tstep]

        if tstep > 0:
            powerL[tstep + 1] = - gammaL * m1 * v1[tstep] * 0.5 * (v1[tstep] + \
                    v1[tstep + 1]) + np.sqrt(2 * kB * gammaL * TL * m1 / dt) * \
                    xiL * 0.5 * (v1[tstep] + v1[tstep + 1])
            powerR[tstep + 1] = - gammaR * m1 * v1[tstep] * 0.5 * (v1[tstep] + \
                    v1[tstep + 1]) + np.sqrt(2 * kB * gammaR * TR * m1 / dt) * \
                    xiR * 0.5 * (v1[tstep] + v1[tstep + 1])
            powerLsq[tstep + 1] = powerL[tstep + 1] * powerL[tstep + 1]
            powerRsq[tstep + 1] = powerR[tstep + 1] * powerR[tstep + 1]

        K1[tstep] = 0.5 * m1 * v1[tstep] ** 2
        K[tstep] = K1[tstep]

        tstep += 1

    return x1, K1, U, K, powerL, powerR, powerLsq, powerRsq 


if __name__ == '__main__':

    SLOTS = int(os.getenv('NSLOTS')) # Get the NSLOTS environment variable provided by the scheduler

    start_time = time.time()

    tBegin = 0.
    tEnd = 1000
    dt = 0.001
    tArray = np.arange(tBegin, tEnd, dt)
    tsize = len(tArray)
    x_t = np.zeros(tsize)
    Utraj = np.zeros(tsize)
    Ktraj = np.zeros(tsize)
    K1traj = np.zeros(tsize)
    powerLtraj = np.zeros(tsize)
    powerRtraj = np.zeros(tsize)
    powerLsqtraj = np.zeros(tsize)
    powerRsqtraj = np.zeros(tsize)

    m1 = float(sys.argv[1])
    mL = 32
    mR = 32
    epsilon = 0.3
    sigma = 3.5
    AL = 1 * 1e6 #  1eV = 23kcal/mol
    alphaL = 5e0
    AR = 1 * 1e6  #  1eV = 23kcal/mol
    alphaR = 5e0
    kB = 0.001987
    TL = 270
    TR = 200
    gammaL = 1.076581 / (m1 / 2)
    gammaR = 1.076581 / (m1 / 2)

    traj = 500


    p = Pool(processes=SLOTS)# pass the number of core to the Pool so that I know how many cores I can use.

    for x, K1, U, K, powerL, powerR, powerLsq, powerRsq in p.map(single_traj, range(traj)):
        x_t += x
        K1traj += K1
        Utraj += U
        Ktraj += K
        powerLtraj += powerL
        powerRtraj += powerR
        powerLsqtraj += powerLsq
        powerRsqtraj += powerRsq

    p.close()
    p.join()


    Utraj /= traj
    Ktraj /= traj
    K1traj /= traj
    powerLtraj /= traj
    powerRtraj /= traj
    powerLsqtraj /= traj
    powerRsqtraj /= traj

    ##----------
    NN = tsize / 2
    # Ksteady = Ktraj[30000:]
    # Kaver = np.sum(Ksteady)/len(Ksteady)
    # print Kaver
    K1steady = K1traj[NN:]
    K1steady = np.mean(K1steady.reshape(-1, 500), axis=1)  # a very neat answer from StackOverflow
    T1aver = np.mean(K1steady) * 2 / kB
    T1std = np.std(K1steady) * 2 / kB
    print 'T1 = ', T1aver, T1std

    PsteadyL = powerLtraj[NN:]
    PsteadyL = np.mean(PsteadyL.reshape(-1, 500), axis=1)  # a very neat answer from StackOverflow
    # to average over a length of numbers
    PsqsteadyL = powerLsqtraj[NN:]
    JLaver = np.mean(PsteadyL)
    JLstd = np.std(PsteadyL)
    JLstd_true = np.sqrt(np.mean(PsqsteadyL) - JLaver**2)
    print 'heatL = ', JLaver, JLstd, JLstd_true

    PsteadyR = powerRtraj[NN:]
    PsteadyR = np.mean(PsteadyR.reshape(-1, 500), axis=1)
    PsqsteadyR = powerRsqtraj[NN:]
    JRaver = np.mean(PsteadyR)
    JRstd = np.std(PsteadyR)
    JRstd_true = np.sqrt(np.mean(PsqsteadyR) - JRaver**2)
    print 'heatR = ', JRaver, JRstd, JRstd_true

    run_time = time.time() - start_time
    print 'run time is: ', run_time / 60.

    ##-----------write-data-out---------
    filename = 'center-' + str(m1) + time.strftime('-%m-%d-%H%M%S.txt')
    with open(filename, "w") as f:
        f.write("time spent in minutes: %f\n" %(run_time/60))
        f.write("AL = %f, alphaL = %f\n" %(AL, alphaL))
        f.write("mass = %f\n" %(m1))
        f.write("trajectory number: %d\n" %(traj))
        f.write("time_step: %f\n" %(dt))
        f.write("number of steps: %d\n" %(tsize/2))
        f.write("TL = %d, TR = %d\n" %(TL, TR))
        f.write("T1 = %f, T1std = %f\n" %(T1aver, T1std))
        f.write("JL = %f, STDJL = %f, STDJL_r = %f\n" %(JLaver, JLstd, JLstd_true))
        f.write("JR = %f, STDJR = %f, STDJR_r = %f\n" %(JRaver, JRstd, JRstd_true))

# filename2 = 'heatflux-singleatom-' + str(str(m1)) + time.strftime('-%m-%d-%H%M%S.txt')
# np.savetxt(filename2, np.c_[PsteadyL, PsteadyR])

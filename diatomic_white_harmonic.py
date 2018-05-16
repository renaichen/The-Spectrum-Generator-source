import sys
import os

import numpy as np
import time
import random
from multiprocessing import Pool

#------------------------------------------------------

def diatomic_traj(n):

    x1 = np.zeros(tsize)
    x2 = np.zeros(tsize)
    xL = np.zeros(tsize)
    xR = np.zeros(tsize)

    v1 = np.zeros(tsize)
    v2 = np.zeros(tsize)

    U = np.zeros(tsize)
    K = np.zeros(tsize)
    K1 = np.zeros(tsize)
    K2 = np.zeros(tsize)
    UintL = np.zeros(tsize)
    UintR = np.zeros(tsize)

    fint = np.zeros(tsize)
    fLt = np.zeros(tsize)
    fRt = np.zeros(tsize)
    powerL = np.zeros(tsize)
    powerR = np.zeros(tsize)
    power12 = np.zeros(tsize)
    powerLsq = np.zeros(tsize)
    powerRsq = np.zeros(tsize)
    power12sq = np.zeros(tsize)

    xL[0] = 46
    x1[0] = 49.
    x2[0] = 51.
    xR[0] = 54.
    halfd = 0.5 * (x2[0] - x1[0])

    v1[0] = 0.
    v2[0] = 0.0

    f1new = 0.0
    f2new = 0.0
    fL = 0.0
    fR = 0.0

    tstep = 0
    while tstep < (tsize-1):
#
        f1old = f1new
        f2old = f2new
        fLold = fL
        fRold = fR

        xiL = random.gauss(0, 1)
        xiR = random.gauss(0, 1)

# #-----------EOM integrator: using the velocity verlet algorithm-----------------------
        x1[tstep+1] = x1[tstep] + v1[tstep]*dt + (0.5/m1)*f1old*dt**2
        x2[tstep+1] = x2[tstep] + v2[tstep]*dt + (0.5/m2)*f2old*dt**2
        xL[tstep+1] = xL[0]
        xR[tstep+1] = xR[0]
        f1new = k12*(x2[tstep+1]-x1[tstep+1]-x012)
        f2new = -f1new
        fint[tstep + 1] = f1new

        # fL = - AL * alphaL * np.exp(-alphaL * (x1[tstep + 1] - xL[tstep + 1]))
        # fR = AR * alphaR * np.exp(-alphaR * (xR[tstep + 1] - x2[tstep + 1]))

        fL = k1l * (x1[tstep+1] - xL[tstep+1] - x1l)
        fR = - k2r * (xR[tstep+1] - x2[tstep+1] - x2r)

        fLt[tstep] = fL
        fRt[tstep] = fR

        if x1[tstep+1] < xL[tstep+1] or x2[tstep+1] > xR[tstep+1] or \
                x1[tstep+1] > x2[tstep+1]:
            f1 = open('wrong-dia-' + str(omega1) + time.strftime('-%m-%d-%H%M%S.txt'), 'w')
            print >> f1, 'error: position disorder, exiting...', omega1, \
                xL[tstep+1], x1[tstep+1], x2[tstep+1], xR[tstep+1] 
            f1.close()
            break

        f1new -= fL
        f2new -= fR
#
        v1[tstep+1] = v1[tstep] + 0.5*((f1old+f1new)/m1) * dt - \
                        gammaL * v1[tstep] * dt + np.sqrt(2 * kB * gammaL * TL * dt / m1) * xiL
        v2[tstep+1] = v2[tstep] + 0.5*((f2old+f2new)/m2) * dt - \
                        gammaR * v2[tstep] * dt + np.sqrt(2 * kB * gammaR * TR * dt / m2) * xiR
# #----------------------------------------------------------------------------------
        U[tstep] = 0.5*k12*(x2[tstep]-x1[tstep]-x012)**2
        UintL[tstep] = AL * np.exp(-alphaL * (x1[tstep] - xL[tstep]))
        UintR[tstep] = AR * np.exp(-alphaR * (xR[tstep] - x2[tstep]))
        U[tstep] += UintL[tstep]
        U[tstep] += UintR[tstep]

        if tstep > 0:
            powerL[tstep + 1] = - gammaL * m1 * v1[tstep] * 0.5 * (v1[tstep] + \
                    v1[tstep + 1]) + np.sqrt(2 * kB * gammaL * TL * m1 / dt) * \
                    xiL * 0.5 * (v1[tstep] + v1[tstep + 1])
            powerR[tstep + 1] = - gammaR * m2 * v2[tstep] * 0.5 * (v2[tstep] + \
                    v2[tstep + 1]) + np.sqrt(2 * kB * gammaR * TR * m2 / dt) * \
                    xiR * 0.5 * (v2[tstep] + v2[tstep + 1])
            power12[tstep + 1] = - 0.5 * fint[tstep + 1] * (v2[tstep + 1] + v1[tstep + 1])
            powerLsq[tstep + 1] = powerL[tstep + 1] * powerL[tstep + 1]
            powerRsq[tstep + 1] = powerR[tstep + 1] * powerR[tstep + 1]
            power12sq[tstep + 1] = power12[tstep + 1] * power12[tstep + 1]

        K1[tstep] = 0.5 * m1 * v1[tstep] ** 2
        K2[tstep] = 0.5 * m2 * v2[tstep] ** 2
        K[tstep] = K1[tstep] + K2[tstep]

        tstep += 1

    return x1, x2, K1, K2,  U, K, powerL, power12, powerR, powerLsq, power12sq, powerRsq 

if __name__ == '__main__':

    SLOTS = int(os.getenv('NSLOTS')) # Get the NSLOTS environment variable provided by the scheduler

    start_time = time.time()

    tBegin = 0.
    tEnd = 1000
    dt = 0.001
    tArray = np.arange(tBegin, tEnd, dt)
    tsize = len(tArray)
    x1_t = np.zeros(tsize)
    x2_t = np.zeros(tsize)
    Utraj = np.zeros(tsize)
    Ktraj = np.zeros(tsize)
    K1traj = np.zeros(tsize)
    K2traj = np.zeros(tsize)
    powerLtraj = np.zeros(tsize)
    powerRtraj = np.zeros(tsize)
    power12traj = np.zeros(tsize)
    powerLsqtraj = np.zeros(tsize)
    powerRsqtraj = np.zeros(tsize)
    power12sqtraj = np.zeros(tsize)

    # m1 = float(sys.argv[1]) / 2.
    # m2 = float(sys.argv[1]) / 2.
    m1 = 24 / 2.
    m2 = 24 / 2.
    mu = m1 * m2 / (m1 + m2)
    mL = 32
    mR = 32
    omegaD = 1.
    omegaL = 1.
    omegaR = 1.
    omega1 = float(sys.argv[1])
    k12 = mu * omega1**2
    k1l = m1 * omegaL**2
    k2r = m2 * omegaR**2
    # x10 = 50
    x012 = 2.0
    x1l = 3.0
    x2r = 3.0
    traj = 0
    epsilon = 0.3
    sigma = 3.5
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

    traj = 40

    p = Pool(processes=SLOTS)# pass the number of core to the Pool so that I know how many cores I can use.

    for x1, x2, K1, K2,  U, K, powerL, powerR, power12,  powerLsq, power12sq, powerRsq in p.map(diatomic_traj, range(traj)):
        x1_t += x1
        x2_t += x2
        K1traj += K1
        K2traj += K2
        Utraj += U
        Ktraj += K
        powerLtraj += powerL
        power12traj += power12
        powerRtraj += powerR
        powerLsqtraj += powerLsq
        power12sqtraj += power12sq
        powerRsqtraj += powerRsq

    p.close()
    p.join()

    Utraj /= traj
    Ktraj /= traj
    K1traj /= traj
    K2traj /= traj
    powerLtraj /= traj
    power12traj /= traj
    powerRtraj /= traj
    powerLsqtraj /= traj
    power12sqtraj /= traj
    powerRsqtraj /= traj
    ##----------
    run_time = time.time() - start_time
    print 'run time is: ', run_time / 60.
    NN = tsize / 2
    # Ksteady = Ktraj[30000:]
    # Kaver = np.sum(Ksteady)/len(Ksteady)
    # print Kaver
    K1steady = K1traj[NN:]
    K1steady = np.mean(K1steady.reshape(-1, 500), axis=1)  # a very neat answer from StackOverflow
    T1aver = np.mean(K1steady) * 2 / kB
    T1std = np.std(K1steady) * 2 / kB
    print 'T1 = ', T1aver, T1std
    K2steady = K2traj[NN:]
    K2steady = np.mean(K2steady.reshape(-1, 500), axis=1)  # a very neat answer from StackOverflow
    T2aver = np.mean(K2steady) * 2 / kB
    T2std = np.std(K2steady) * 2 / kB
    print 'T2 = ', T2aver, T2std


    PsteadyL = powerLtraj[NN:]
    PsteadyL = np.mean(PsteadyL.reshape(-1, 500), axis=1)  # a very neat answer from StackOverflow
    # to average over a length of numbers
    PsqsteadyL = powerLsqtraj[NN:]
    JLaver = np.mean(PsteadyL)
    JLstd = np.std(PsteadyL)
    JLstd_true = np.sqrt(np.mean(PsqsteadyL) - JLaver**2)
    print 'heatL = ', JLaver, JLstd, JLstd_true
    Psteady12 = power12traj[NN:]
    Psteady12 = np.mean(Psteady12.reshape(-1, 500), axis=1)
    Psqsteady12 = power12sqtraj[NN:]
    J12aver = np.mean(Psteady12)
    J12std = np.std(Psteady12)
    J12std_true = np.sqrt(np.mean(Psqsteady12) - J12aver**2)
    print 'heat12 = ', J12aver, J12std, J12std_true
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
filename = 'diatomic-' + str(omega1) + time.strftime('-%m-%d-%H%M%S.txt')
with open(filename, "w") as f:
    f.write("time spent in minutes: %f\n" %(run_time/60))
    f.write("AL = %f, alphaL = %f\n" %(AL, alphaL))
    f.write("mass = %f\n" %(m1))
    f.write("omega_r = %f\n" %(omega1))
    f.write("equilibrium length (bond length): %f\n" %(x012))
    f.write("trajectory number: %d\n" %(traj))
    f.write("time_step: %f\n" %(dt))
    f.write("number of steps: %d\n" %(tsize/2))
    f.write("TL = %d, TR = %d\n" %(TL, TR))
    f.write("T1 = %f, T1std = %f\n" %(T1aver, T1std))
    f.write("T2 = %f, T2std = %f\n" %(T2aver, T2std))
    f.write("JL = %f, STDJL = %f, STDJL_r = %f\n" %(JLaver, JLstd, JLstd_true))
    f.write("J12 = %f, STDJ12 = %f, STDJ12_r = %f\n" %(J12aver, J12std, J12std_true))
    f.write("JR = %f, STDJR = %f, STDJR_r = %f\n" %(JRaver, JRstd, JRstd_true))

# filename2 = 'heatflux-diatomic-' + str(m1) + time.strftime('-%m-%d-%H%M%S.txt')
# np.savetxt(filename2, np.c_[PsteadyL, Psteady12, PsteadyR])


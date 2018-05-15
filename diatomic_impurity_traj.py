import sys
import os

import numpy as np
import time
import Debye_spectrum_3 as ds
from multiprocessing import Pool

#------------------------------------------------------

def diatomic_singletraj(n):
    sp_objL = ds.Generator(n=nL,
                           mass=massL,
                           omegaD=omegaDL,
                           temperature=temperatureL,
                           dt=dt1L,
                           t_num=t_numL,
                           Ntraj=Ntraj1L,
                           )
    rand_arrayL = sp_objL.give_me_random_series(dt)

    sp_objR = ds.Generator(n=nR,
                           mass=massR,
                           omegaD=omegaDR,
                           temperature=temperatureR,
                           dt=dt1R,
                           t_num=t_numR,
                           Ntraj=Ntraj1R,
                           )
    rand_arrayR = sp_objR.give_me_random_series(dt)

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
    term1 = np.zeros(tsize)
    term2 = np.zeros(tsize)

    damperL = np.zeros(tsize)
    damperR = np.zeros(tsize)

    xL[0] = 46
    x1[0] = 49.
    x2[0] = 51
    xR[0] = 54
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
        damperL[tstep] = sp_objL.damp_getter(fL)
        damperR[tstep] = sp_objR.damp_getter(fR)

# #-----------EOM integrator: using the velocity verlet algorithm-----------------------
        x1[tstep+1] = x1[tstep] + v1[tstep]*dt + (0.5/m1)*f1old*dt**2
        x2[tstep+1] = x2[tstep] + v2[tstep]*dt + (0.5/m2)*f2old*dt**2
        xL[tstep+1] = xL[0] + damperL[tstep] + rand_arrayL[tstep]
        xR[tstep+1] = xR[0] + damperR[tstep] + rand_arrayR[tstep]
        f1new = k12*(x2[tstep+1]-x1[tstep+1]-x012)
        f2new = -f1new
        fint[tstep + 1] = f1new

        fL = - AL * alphaL * np.exp(-alphaL * (x1[tstep + 1] - xL[tstep + 1]))
        fR = AR * alphaR * np.exp(-alphaR * (xR[tstep + 1] - x2[tstep + 1]))
        fLt[tstep] = fL
        fRt[tstep] = fR

        # f2old = 48 * epsilon * sigma**12 / (x1[tstep+1]-x2[tstep+1])**13 \
        #             - 24 * epsilon * sigma**6/(x1[tstep+1]-x2[tstep+1])**7
        # f2old = 10/(x1[tstep+1]-x2[tstep+1])
        f1new -= fL
        f2new -= fR
#
        v1[tstep+1] = v1[tstep] + 0.5*((f1old+f1new)/m1) * dt
        v2[tstep+1] = v2[tstep] + 0.5*((f2old+f2new)/m2) * dt
# #----------------------------------------------------------------------------------
        U[tstep] = 0.5*k12*(x2[tstep]-x1[tstep]-x012)**2
        UintL[tstep] = AL * np.exp(-alphaL * (x1[tstep] - xL[tstep]))
        UintR[tstep] = AR * np.exp(-alphaR * (xR[tstep] - x2[tstep]))
        U[tstep] += UintL[tstep]
        U[tstep] += UintR[tstep]

        # powerL[tstep + 1] = fint[tstep + 1] * v1[tstep + 1] + (UintL[tstep + 1] - UintL[tstep]) / dt
        # powerR[tstep + 1] = -fint[tstep + 1] * v2[tstep + 1] + (UintR[tstep + 1] - UintR[tstep]) / dt
        # # powerL[tstep] = f1old * v1[tstep] + 0.5 * (UintL[tstep] - UintL[tstep - 1])/dt
        # # powerR[tstep] = f1old * v1[tstep] + 0.5 * (UintR[tstep] - UintR[tstep - 1])/dt
        # power12[tstep + 1] = 0.5 * fint[tstep + 1] * (v2[tstep + 1] - v1[tstep + 1])
        # term1[tstep + 1] = fint[tstep + 1] * v1[tstep + 1]
        # term2[tstep + 1] = (UintL[tstep + 1] - UintL[tstep]) / dt

        if tstep > 0:
            powerL[tstep] = - 0.5 * fL * ((xL[tstep + 1] - xL[tstep]) / dt + v1[tstep])
            powerR[tstep] = 0.5 * fR * ((xR[tstep + 1] - xR[tstep]) / dt + v2[tstep])
            power12[tstep] = - 0.5 * fint[tstep] * (v2[tstep] + v1[tstep])
            powerLsq[tstep] = powerL[tstep] * powerL[tstep]
            powerRsq[tstep] = powerR[tstep] * powerR[tstep]
            power12sq[tstep] = power12[tstep] * power12[tstep]
            term1[tstep] = fint[tstep] * v1[tstep]

        # U[tstep] += 4 * epsilon * sigma**12 / (x1[tstep]-x2[tstep])**12\
        #                 - 4 * epsilon * sigma**6/(x1[tstep]-x2[tstep])**6
        K1[tstep] = 0.5 * m1 * v1[tstep] ** 2
        K2[tstep] = 0.5 * m2 * v2[tstep] ** 2
        K[tstep] = K1[tstep] + K2[tstep]

        tstep += 1
  #      print "time steps: ", tstep

    return x1, x2, U, K1, K2, K, powerL, power12, powerR, powerLsq, power12sq, powerRsq


if __name__=='__main__':

    SLOTS = int(os.getenv('NSLOTS')) # Get the NSLOTS environment variable provided by the scheduler

    start_time = time.time()
    
    tBegin = 0.0
    tEnd = 100
    dt = 0.001
    kB = 0.00198
    Tpoint = 200000

    nL = 8
    omegaDL = 2.
    massL = 32.
    t_numL = Tpoint
    dt1L = dt
    Ntraj1L = 1
    temperatureL = 450.0

    nR = 8
    omegaDR = 2.
    massR = 32.
    t_numR = Tpoint
    dt1R = dt
    Ntraj1R = 1
    temperatureR = 250.0

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
    damper_trajL = np.zeros(tsize)
    damper_trajR = np.zeros(tsize)

    m1 = float(sys.argv[1]) / 2.
    m2 = float(sys.argv[1]) / 2.
    mu = m1 * m2 / (m1 + m2)
    omega1 = 1.
    k12 = mu * omega1**2
    # x10 = 50
    x012 = 2.0
    epsilon = 0.3
    sigma = 3.5
    AL = 1 * 1e6  #  1eV = 23kcal/mol
    alphaL = 5e0
    AR = 1 * 1e6  #  1eV = 23kcal/mol
    alphaR = 5e0

    traj = 45

    p = Pool(processes=SLOTS)# pass the number of core to the Pool so that I know how many cores I can use.

    for x1, x2, U, K1, K2, K, powerL, power12, powerR, powerLsq, power12sq, powerRsq in p.map(diatomic_singletraj, range(traj)):
        x1_t += x1
        x2_t += x2
        Utraj += U
        K1traj += K1
        K2traj += K2
        Ktraj += K
        powerLtraj += powerL
        power12traj += power12
        powerRtraj += powerR
        powerLsqtraj += powerLsq
        power12sqtraj += power12sq
        powerRsqtraj += powerRsq

    p.close()
    p.join()


    x1_t /= traj
    x2_t /= traj
    Utraj /= traj
    K1traj /= traj
    K2traj /= traj
    Ktraj /= traj
    powerLtraj /= traj
    powerRtraj /= traj
    power12traj /= traj
    powerLsqtraj /= traj
    powerRsqtraj /= traj
    power12sqtraj /= traj

##----------
    NN = tsize / 2
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

    ##-----------write-data-out---------
    #filename = time.strftime('diatomic-%m-%d-%H%M%S.txt')
    filename = 'diatomic-' + str(m1) + time.strftime('-%m-%d-%H%M%S.txt')
    with open(filename, "w") as f:
        f.write("time spent in minutes: %f\n" %(run_time/60))
        f.write("AL = %f, alphaL = %f\n" %(AL, alphaL))
        f.write("mass = %f\n" %(m1))
        f.write("equilibrium length (bond length): %f\n" %(x012))
        # f.write("initial postitions: xL = %f, x1 = %f, x2 = %f, xR = %f\n" %(xL[0], x1[0], x2[0], xR[0]))
        f.write("omegaD = %f, omega = %f\n" %(omegaDL, omega1))
        f.write("trajectory number: %d\n" %(traj))
        f.write("time_step: %f\n" %(dt))
        f.write("number of steps: %d\n" %(tsize/2))
        f.write("TL = %d, TR = %d\n" %(temperatureL, temperatureR))
        f.write("T1 = %f, T1std = %f\n" %(T1aver, T1std))
        f.write("T2 = %f, T2std = %f\n" %(T2aver, T2std))
        f.write("JL = %f, STDJL = %f, STDJL_r = %f\n" %(JLaver, JLstd, JLstd_true))
        f.write("J12 = %f, STDJ12 = %f, STDJ12_r = %f\n" %(J12aver, J12std, J12std_true))
        f.write("JR = %f, STDJR = %f, STDJR_r = %f\n" %(JRaver, JRstd, JRstd_true))

    #filename2 = time.strftime('heatflux-diatomic-%m-%d-%H%M%S.txt')
    filename2 = 'heatflux-diatomic-' + str(m1) + time.strftime('-%m-%d-%H%M%S.txt')
    np.savetxt(filename2, np.c_[PsteadyL, Psteady12, PsteadyR])


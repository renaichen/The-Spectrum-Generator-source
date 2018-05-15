import sys
import os

import numpy as np
import Debye_spectrum_3 as ds
import time
from multiprocessing import Pool

#------------------------------------------------------

def single_singletraj(n):
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

    damperL = np.zeros(tsize)
    damperR = np.zeros(tsize)

    xL[0] = 46
    x1[0] = 50.
    xR[0] = 54
    halfd = 1.

    v1[0] = 0.0

    f1new = 0.0
    fL = 0.0
    fR = 0.0

    tstep = 0
    while tstep < (tsize-1):
#
        f1old = f1new
        damperL[tstep] = sp_objL.damp_getter(fL)
        damperR[tstep] = sp_objR.damp_getter(fR)

# #-----------EOM integrator: using the velocity verlet algorithm-----------------------
        x1[tstep+1] = x1[tstep] + v1[tstep]*dt + (0.5/m1)*f1old*dt**2
        xL[tstep+1] = xL[0] + damperL[tstep] + rand_arrayL[tstep]
        xR[tstep+1] = xR[0] + damperR[tstep] + rand_arrayR[tstep]
        # xL[tstep + 1] = xL[0] + rand_arrayL[tstep]
        # xR[tstep + 1] = xR[0] + rand_arrayR[tstep]

        fL = - AL * alphaL * np.exp(-alphaL * (x1[tstep + 1] - halfd - xL[tstep + 1]))
        fR = AR * alphaR * np.exp(-alphaR * (xR[tstep + 1] - halfd - x1[tstep + 1]))

        if x1[tstep+1] < xL[tstep+1] or x2[tstep+1] > xR[tstep+1] or \
                x1[tstep+1] > x2[tstep+1]:
            f1 = open('wrong-dia-' + str(m1) + time.strftime('-%m-%d-%H%M%S.txt'), 'w')
            print >> f1, 'error: position disorder, exiting...', omega1, \
                xL[tstep+1], x1[tstep+1], x2[tstep+1], xR[tstep+1] 
            f1.close()
            break

        fLt[tstep] = fL
        fRt[tstep] = fR

        f1new = - (fL + fR)
        # if tstep < 1000:
        #     print f1new, x1[tstep]
#
        v1[tstep+1] = v1[tstep] + 0.5*((f1old+f1new)/m1) * dt
# #----------------------------------------------------------------------------------
        UintL[tstep] = AL * np.exp(-alphaL * (x1[tstep] - halfd - xL[tstep]))
        UintR[tstep] = AR * np.exp(-alphaR * (xR[tstep] - halfd - x1[tstep]))
        U[tstep] += UintL[tstep]
        U[tstep] += UintR[tstep]

        if tstep > 0:
            powerL[tstep] = - 0.5 * fL * ((xL[tstep + 1] - xL[tstep]) / dt + v1[tstep])
            powerR[tstep] = 0.5 * fR * ((xR[tstep + 1] - xR[tstep]) / dt + v1[tstep])
            powerLsq[tstep] = powerL[tstep] * powerL[tstep]
            powerRsq[tstep] = powerR[tstep] * powerR[tstep]

        K1[tstep] = 0.5 * m1 * v1[tstep] ** 2
        K[tstep] = K1[tstep]

        tstep += 1

    return x1, K1, U, K, powerL, powerR, powerLsq, powerRsq 



if __name__ == '__main__':

    SLOTS = int(os.getenv('NSLOTS')) # Get the NSLOTS environment variable provided by the scheduler

    start_time = time.time()

    tBegin = 0.0
    tEnd = 250.
    dt = 0.005
    kB = 0.00198
    Tpoint = 100000

    nL = 8
    omegaDL = 2
    massL = 32.
    t_numL = Tpoint
    dt1L = dt
    Ntraj1L = 1
    temperatureL = 50.0

    nR = 8
    omegaDR = 2
    massR = 32.
    t_numR = Tpoint
    dt1R = dt
    Ntraj1R = 1
    temperatureR = 50.0

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
    damper_trajL = np.zeros(tsize)
    damper_trajR = np.zeros(tsize)

    m1 = float(sys.argv[1])
    traj = 5000
    AL = 1 * 1e6 #  1eV = 23kcal/mol
    alphaL = 5e0
    AR = 1 * 1e6  #  1eV = 23kcal/mol
    alphaR = 5e0

    p = Pool(processes=SLOTS)# pass the number of core to the Pool so that I know how many cores I can use.

    for x, K1, U, K, powerL, powerR, powerLsq, powerRsq in p.map(single_singletraj, range(traj)):
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
    #filename = time.strftime('center-atom-%m-%d-%H%M%S.txt')
    filename = 'center-' + str(m1) + time.strftime('-%m-%d-%H%M%S.txt')
    with open(filename, "w") as f:
        f.write("time spent in minutes: %f\n" %(run_time/60))
        f.write("AL = %f, alphaL = %f\n" %(AL, alphaL))
        f.write("mass = %f\n" %(m1))
        # f.write("initial postitions: xL = %f, x1 = %f, xR = %f\n" %(xL[0], x1[0], xR[0]))
        f.write("trajectory number: %d\n" %(traj))
        f.write("time_step: %f\n" %(dt))
        f.write("number of steps: %d\n" %(tsize/2))
        f.write("TL = %d, TR = %d\n" %(temperatureL, temperatureR))
        f.write("T1 = %f, T1std = %f\n" %(T1aver, T1std))
        # f.write("T2 = %f\n" %(T2aver))
        f.write("JL = %f, STDJL = %f, STDJL_r = %f\n" %(JLaver, JLstd, JLstd_true))
        f.write("JR = %f, STDJR = %f, STDJR_r = %f\n" %(JRaver, JRstd, JRstd_true))

    select_out = tsize / 200
    new_time = tArray[0:-1:select_out]
    new_K1traj = K1traj[0:-1:select_out]

    filename2 = 'kinetic-singleatom-' + str(str(m1)) + time.strftime('-%m-%d-%H%M%S.txt')
    np.savetxt(filename2, np.c_[new_time, new_K1traj, new_K1traj * 2 / kB])

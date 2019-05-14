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

    K1 = 0.0

    powerL = 0.0
    powerR = 0.0

    xLeq = 46.
    x1new = 50.
    xReq = 54.
    xLnew = xLeq
    xRnew = xReq

    v1new = 0.0

    f1new = 0.0
    fL = 0.0
    fR = 0.0

    tstep = 0
    while tstep < (tsize-1):
#
        f1old = f1new
        fLold = fL
        fRold = fR
        x1old = x1new
        xLold = xLnew
        xRold = xRnew
        v1old = v1new
        damperL = sp_objL.damp_getter(fL)
        damperR = sp_objR.damp_getter(fR)

# #-----------EOM integrator: using the velocity verlet algorithm-----------------------
        x1new = x1old + v1old*dt + (0.5/m1)*f1old*dt**2
        xLnew = xLeq + damperL + rand_arrayL[tstep]
        xRnew = xReq + damperR + rand_arrayR[tstep]

        fL = - AL * alphaL * np.exp(-alphaL * (x1new - halfd - xLnew))
        fR = AR * alphaR * np.exp(-alphaR * (xRnew - -halfd - x1new))

        # fL = k1l * (x1new - xL - x1l)
        # fR = - k2r * (xR - x2new - x2r)

        # fL = - AL * alphaL * np.exp(-alphaL * (x1[tstep + 1] - halfd - xL[tstep + 1]))
        # fR = AR * alphaR * np.exp(-alphaR * (xR[tstep + 1] - halfd - x1[tstep + 1]))

        if x1new < xLnew or x1new > xRnew:
            f1 = open('wrong-single-' + time.strftime('-%m-%d-%H%M%S.txt'), 'w')
            print >> f1, 'error: position disorder, exiting...', omega1, \
                xLnew, x1new, xRnew 
            f1.close()
            break

        f1new = - (fL + fR)
#
        v1new = v1old + 0.5*((f1old+f1new)/m1) * dt
# #----------------------------------------------------------------------------------
        # UintL[tstep] = AL * np.exp(-alphaL * (x1[tstep] - halfd - xL[tstep]))
        # UintR[tstep] = AR * np.exp(-alphaR * (xR[tstep] - halfd - x1[tstep]))
        # U[tstep] += UintL[tstep]
        # U[tstep] += UintR[tstep]

        if tstep > (tsize / 2 - 1):
            powerL += - 0.5 * fL * ((xLnew - xLold) / dt + v1old)
            powerR += 0.5 * fR * ((xRnew - xRold) / dt + v1old)
            K1 += 0.5 * m1 * v1old ** 2


        tstep += 1

    return K1, powerL, powerR 


if __name__ == '__main__':

    # SLOTS = int(os.getenv('NSLOTS')) # Get the NSLOTS environment variable provided by the scheduler
    SLOTS = int(os.getenv('SLURM_NPROCS')) # For SLURM

    start_time = time.time()

    traj = 20

    tBegin = 0.0
    tEnd = 100.
    dt = 0.001
    kB = 0.00198
    Tpoint = 200000

    nL = 8
    omegaDL = 2
    massL = 32.
    t_numL = Tpoint
    dt1L = dt
    Ntraj1L = 1
    temperatureL = 270.0

    nR = 8
    omegaDR = 2
    massR = 32.
    t_numR = Tpoint
    dt1R = dt
    Ntraj1R = 1
    temperatureR = 200.0

    tArray = np.arange(tBegin, tEnd, dt)
    tsize = len(tArray)

    halftsize = tsize / 2
    Ktraj = np.zeros(traj)
    K1traj = np.zeros(traj)
    K2traj = np.zeros(traj)
    powerLtraj = np.zeros(traj)
    powerRtraj = np.zeros(traj)
    power12traj = np.zeros(traj)

    m1 = float(sys.argv[1])
    AL = 1 * 1e6 #  1eV = 23kcal/mol
    alphaL = 5e0
    AR = 1 * 1e6  #  1eV = 23kcal/mol
    alphaR = 5e0
    halfd = 1.0
    
    p = Pool(processes=SLOTS)# pass the number of core to the Pool so that I know how many cores I can use.

    i = 0
    for K1, powerL, powerR in p.map(single_singletraj, range(traj)):
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
    #filename = time.strftime('center-atom-%m-%d-%H%M%S.txt')
    filename = 'center-' + str(m1) + time.strftime('-%m-%d-%H%M%S.txt')
    with open(filename, "w") as f:
        f.write("time spent in minutes: %f\n" %(run_time/60))
        # f.write("AL = %f, alphaL = %f\n" %(AL, alphaL))
        f.write("mass = %f\n" %(m1))
        # f.write("equilibrium length (bond length): %f\n" %(x012))
        f.write("trajectory number: %d\n" %(traj))
        f.write("time_step: %f\n" %(dt))
        f.write("number of steps: %d\n" %(tsize/2))
        f.write("TL = %d, TR = %d\n" %(temperatureL, temperatureR))
        f.write("T1 = %f, T1std = %f\n" %(T1aver, T1std))
        f.write("JL = %f, STDJL = %f\n" %(JLaver, JLstd))
        f.write("JR = %f, STDJR = %f\n" %(JRaver, JRstd))

    # filename2 = 'kinetic-singleatom-' + str(str(m1)) + time.strftime('-%m-%d-%H%M%S.txt')
    # np.savetxt(filename2, np.c_[new_time, new_K1traj, new_K1traj * 2 / kB])

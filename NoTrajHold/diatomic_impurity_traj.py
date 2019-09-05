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


    K = 0.0
    K1 = 0.0
    K2 = 0.0

    fint = 0.0
    powerL = 0.0
    powerR = 0.0
    power12 = 0.0

    xLeq = 46
    x1new = 49.
    x2new = 51
    xReq = 54
    xLnew = xLeq
    xRnew = xReq

    v1new = 0.0
    v2new = 0.0

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
        x1old = x1new
        x2old = x2new
        xLold = xLnew
        xRold = xRnew
        v1old = v1new
        v2old = v2new
        damperL = sp_objL.damp_getter(fL)
        damperR = sp_objR.damp_getter(fR)

# #-----------EOM integrator: using the velocity verlet algorithm-----------------------
        x1new = x1old + v1old*dt + (0.5/m1)*f1old*dt**2
        x2new = x2old + v2old*dt + (0.5/m2)*f2old*dt**2
        xLnew = xLeq + damperL + rand_arrayL[tstep]
        xRnew = xReq + damperR + rand_arrayR[tstep]
        f1new = k12*(x2new-x1new-x012)
        f2new = -f1new
        fint = f1new

        fL = - AL * alphaL * np.exp(-alphaL * (x1new - xLnew))
        fR = AR * alphaR * np.exp(-alphaR * (xRnew - x2new))

        # fL = k1l * (x1new - xL - x1l)
        # fR = - k2r * (xR - x2new - x2r)

        if x1new < xLnew or x2new > xRnew or \
                x1new > x2new:
            f1 = open('wrong-dia-' + str(omega1) + time.strftime('-%m-%d-%H%M%S.txt'), 'w')
            print >> f1, 'error: position disorder, exiting...', omega1, \
                xLnew, x1new, x2new, xRnew 
            f1.close()
            break

        # f2old = 48 * epsilon * sigma**12 / (x1[tstep+1]-x2[tstep+1])**13 \
        #             - 24 * epsilon * sigma**6/(x1[tstep+1]-x2[tstep+1])**7
        # f2old = 10/(x1[tstep+1]-x2[tstep+1])
        f1new -= fL
        f2new -= fR
#
        v1new = v1old + 0.5*((f1old+f1new)/m1) * dt
        v2new = v2old + 0.5*((f2old+f2new)/m2) * dt
# #----------------------------------------------------------------------------------
        # U[tstep] = 0.5*k12*(x2[tstep]-x1[tstep]-x012)**2
        # UintL[tstep] = AL * np.exp(-alphaL * (x1[tstep] - xL[tstep]))
        # UintR[tstep] = AR * np.exp(-alphaR * (xR[tstep] - x2[tstep]))
        # U[tstep] += UintL[tstep]
        # U[tstep] += UintR[tstep]

        if tstep > (tsize / 2 - 1):
            powerL += - 0.5 * fL * ((xLnew - xLold) / dt + v1old)
            powerR += 0.5 * fR * ((xRnew - xRold) / dt + v2old)
            power12 += - 0.5 * fint * (v2old + v1old)

        # U[tstep] += 4 * epsilon * sigma**12 / (x1[tstep]-x2[tstep])**12\
        #                 - 4 * epsilon * sigma**6/(x1[tstep]-x2[tstep])**6
            K1 += 0.5 * m1 * v1old ** 2
            K2 += 0.5 * m2 * v2old ** 2
            K = K1 + K2

        tstep += 1
  #      print "time steps: ", tstep

    return K1, K2, K, powerL, power12, powerR


if __name__=='__main__':

    # SLOTS = int(os.getenv('NSLOTS')) # Get the NSLOTS environment variable provided by the scheduler
    SLOTS = int(os.getenv('SLURM_NPROCS')) # For SLURM

    start_time = time.time()
    
    traj = 23

    tBegin = 0.0
    tEnd = 10000
    dt = 0.01
    kB = 0.00198
    Tpoint = 2000000

    nL = 8
    omegaDL = 2.
    massL = 32.
    t_numL = Tpoint
    dt1L = dt
    Ntraj1L = 1
    temperatureL = 270.0

    nR = 8
    omegaDR = 2.
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

    m1 = float(sys.argv[1]) / 2.
    m2 = float(sys.argv[1]) / 2.
    mu = m1 * m2 / (m1 + m2)
    omega1 = 1.
    k12 = mu * omega1**2
    # x10 = 50
    x012 = 2.0
    # epsilon = 0.3
    # sigma = 3.5
    AL = 1 * 1e6  #  1eV = 23kcal/mol
    alphaL = 5e0
    AR = 1 * 1e6  #  1eV = 23kcal/mol
    alphaR = 5e0

    p = Pool(processes=SLOTS)# pass the number of core to the Pool so that I know how many cores I can use.

## Two ways (with index i or using enumerate by remember using tuple() sign)
## seem all good for indexing the trajectory number, but at the safe side
## using i++ is chosen here.
    i = 0
    for K1, K2, K, powerL, power12, powerR in p.map(diatomic_singletraj, range(traj)):
        K1traj[i] = K1 / halftsize
        K2traj[i] = K2 / halftsize
        Ktraj[i] = K / halftsize
        powerLtraj[i] = powerL / halftsize
        power12traj[i] = power12 / halftsize
        powerRtraj[i] = powerR / halftsize
        i += 1

    p.close()
    p.join()

##----------

    T1aver = np.mean(K1traj) * 2 / kB
    T1std = np.std(K1traj) * 2 / kB
    print 'T1 = ', T1aver, T1std
    T2aver = np.mean(K2traj) * 2 / kB
    T2std = np.std(K2traj) * 2 / kB
    print 'T2 = ', T2aver, T2std

    JLaver = np.mean(powerLtraj)
    JLstd = np.std(powerLtraj)
    print 'heatL = ', JLaver, JLstd
    J12aver = np.mean(power12traj)
    J12std = np.std(power12traj)
    print 'heat12 = ', J12aver, J12std
    JRaver = np.mean(powerRtraj)
    JRstd = np.std(powerRtraj)
    print 'heatR = ', JRaver, JRstd

    run_time = time.time() - start_time
    print 'run time is: ', run_time / 60.

    ##-----------write-data-out---------
    #filename = time.strftime('diatomic-%m-%d-%H%M%S.txt')
    filename = 'diatomic-' + str(m1) + time.strftime('-%m-%d-%H%M%S.txt')
    with open(filename, "w") as f:
        f.write("time spent in minutes: %f\n" %(run_time/60))
        f.write("AL = %f, alphaL = %f\n" %(AL, alphaL))
        f.write("mass = %f\n" %(m1))
        # f.write("equilibrium length (bond length): %f\n" %(x012))
        # f.write("initial postitions: xL = %f, x1 = %f, x2 = %f, xR = %f\n" %(xL[0], x1[0], x2[0], xR[0]))
        f.write("omegaD = %f, omega_r = %f\n" %(omegaDL, omega1))
        f.write("trajectory number: %d\n" %(traj))
        f.write("time_step: %f\n" %(dt))
        f.write("number of steps: %d\n" %(tsize/2))
        f.write("TL = %d, TR = %d\n" %(temperatureL, temperatureR))
        f.write("T1 = %f, T1std = %f\n" %(T1aver, T1std))
        f.write("T2 = %f, T2std = %f\n" %(T2aver, T2std))
        f.write("JL = %f, STDJL = %f\n" %(JLaver, JLstd))
        f.write("J12 = %f, STDJ12 = %f\n" %(J12aver, J12std))
        f.write("JR = %f, STDJR = %f\n" %(JRaver, JRstd))

    # #filename2 = time.strftime('heatflux-diatomic-%m-%d-%H%M%S.txt')
    # filename2 = 'heatflux-diatomic-' + str(m1) + time.strftime('-%m-%d-%H%M%S.txt')
    # np.savetxt(filename2, np.c_[PsteadyL, Psteady12, PsteadyR])


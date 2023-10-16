
"""
Returns a propagating 2D wave of type:
V(r,t) = sum_i [ Vo_i * cos(w_i*t - k_i*r)]


Wave will be saved as a [N x N x t] array with t time values

-dt determined by requirement to resolve the Nyquist freq (plus some enhancement factor), and total
number of timesteps determined to allow wave to propagate across grid.

-number of grid points determined to resolve the smallest wavelength wave

-Amplitude will be distributed as a falling log spectrum with log-spaced frequency values b/t f0 and f1
-k values determined from linear dispersion relation (constant Vph) as 
    f = Vph*k/2pi = Vph/lambda
    Using a constant value of Vph thus gives this wave a range of wavelengths and k-values.
    Assume a constant propagation direction. 


"""


def plane_wave_2d(gridspan=100, #size of grid in meters
                  f0=10, f1=20, #range of freqs (Hz)
                  nfreqs=2,    #number of freqs
                  Vph=1000,     #phase velocity in m/s
                  angle=45,     #angle relative to x-axis (deg)
                  p0=1,p1=0.1,   #power at f0 and f1
                  snapshot=0):  #Single time instance only (t=0)


    import numpy as np 
    import sys 
    sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
    sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
    import matplotlib.pyplot as plt



    #-------------------------------
    #Define timestep to, at minimum, resolve the Nyquist freq.
    #In practice it needs to be quite a bit better than this. 
    timeEnhance = 8
    fnyq = 2*f1
    dt = 1/fnyq / timeEnhance


    #Define number of total time steps to give wave time to propagate across grid. 
    dt_total = gridspan / Vph  #total time needed for wave to propagate across grid 

    if snapshot: 
        ntsteps = 1
    else:
        ntsteps = int(dt_total / dt)
        if ntsteps < 1: 
            ntsteps = 1

    tvals = np.linspace(0,dt_total,num=ntsteps)  
    #-------------------------------



    #Define frequencies (Hz), k-values, and wavelengths
    fvals = np.logspace(np.log10(f0),np.log10(f1),num=nfreqs)
    wvals = 2*np.pi*fvals

    kmag = 2*np.pi*fvals / Vph
    kx = kmag * np.cos(np.radians(angle))
    ky = kmag * np.sin(np.radians(angle))


    #Determine range of wavelengths 
    lmin = 2*np.pi/np.max(kmag)
    lmax = 2*np.pi/np.min(kmag)

    #-------------------------------


    #Number of grid points needs to be, at minimum, high enough to resolve smallest wavelength.
    #Practically, it needs to be quite a bit larger than this. 
    gridEnhance = 32
    ngridpts = int(gridEnhance*(gridspan / lmin))
    xvals = np.linspace(0,np.floor(gridspan),ngridpts)
    yvals = xvals



    #Define a falling spectrum 
    Vo = np.logspace(np.log10(p0), np.log10(p1), nfreqs)





    #--------------------------------------------------
    #Final wave values in array of size [NxNxt]. To avoid nested for loops, broadcast all grid and time arrays to this size 

    


    #grid values for each time
    xvals2 = np.repeat(xvals[:,np.newaxis],ngridpts,axis=1)
    yvals2 = xvals2.T

    xvals2 = np.repeat(xvals2[:,:,np.newaxis],ntsteps,axis=2)
    yvals2 = np.repeat(yvals2[:,:,np.newaxis],ntsteps,axis=2)

    tvals2 = np.repeat(tvals[np.newaxis,:],ngridpts, axis=0)
    tvals2 = np.repeat(tvals2[np.newaxis,:,:],ngridpts, axis=0)




    V = np.zeros((ngridpts, ngridpts, ntsteps))
    for w in range(nfreqs):
        V = V + Vo[w]*np.cos(wvals[w]*tvals2 - (kx[w]*xvals2 + ky[w]*yvals2))



    return(V, xvals, yvals, tvals, fvals, lmin, lmax)




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

NOTE: returned arrays can become very huge if not careful. Best way to avoid this is to make
the grid only a bit larger than the tip-tip boom length of a hypothetical detector. 
Good value is 2x boom length. Smaller values like 1.2 can lead to phasing.

Can also determine Vphase based on desired wavelength of plane wave 
l = 30    # desired wavelength meters
f = 5000  # freq in Hz
Vph = f*l

TESTING: 
-Phase velocity correct (tested by comparing antenna potentials for wavelength >> boom length)
-Animation shows that the wave is propgating at correct angle.
-Wavelength and k-values verified to be correct from snapshot animation
-Frequency correct from single-ended potentials
"""


def plane_wave_2d(gridspan=100, #size of grid in meters
                  f0=10, f1=20, #range of freqs (Hz)
                  nfreqs=2,    #number of freqs
                  Vph=1000,     #phase velocity in m/s (NOTE: the lower Vph the more time steps will be needed to allow for wave to propagate across grid)
                  angle=45,     #angle relative to x-axis (deg)
                  p0=1,p1=0.1,   #power at f0 and f1
                  timeMultiplier=1, #add more runtime for the wave
                  snapshot=0,  #Single time instance only (t=0)
                  animation=0, #watch animation of the wave?
                  no_noise=0,  #don't add in noise spectrum
                  gridEnhance=16, #additional enhancement of spatial grid
                  timeEnhance=16): #additional temporal enhancement


    import numpy as np 
    #import sys 
    #sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
    #sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
    import matplotlib.pyplot as plt
    from perlin_numpy import generate_perlin_noise_3d



    #-------------------------------
    #Define timestep to, at minimum, resolve the Nyquist freq.
    #In practice it needs to be quite a bit better than this. 
    fnyq = 2*f1
    dt = 1/fnyq / timeEnhance


    #Define number of total time steps to give wave time to propagate across grid. 
    dt_total = timeMultiplier* (gridspan / Vph)  #total time needed for wave to propagate across grid 

    if snapshot: 
        ntsteps = 1
    else:
        ntsteps = int(np.round(dt_total / dt))
        if ntsteps < 1: 
            ntsteps = 1

    #want ntsteps to be an even number so that the noise array can be constructed with same size 
    if np.mod(ntsteps, 2):
        ntsteps += 1


    tvals = np.linspace(0,dt_total,num=ntsteps)  
    #-------------------------------



    #Define frequencies (Hz), k-values, and wavelengths
    fvals = np.logspace(np.log10(f0),np.log10(f1),num=nfreqs)
    wvals = 2*np.pi*fvals

    kmag = wvals / Vph   #rad/m
    kx = kmag * np.cos(np.radians(angle))
    ky = kmag * np.sin(np.radians(angle))


    #Determine range of wavelengths 
    lmin = 2*np.pi/np.max(kmag)
    lmax = 2*np.pi/np.min(kmag)

    #-------------------------------


    #Number of grid points needs to be, at minimum, high enough to resolve smallest wavelength.
    #Practically, it needs to be quite a bit larger than this. 
    ngridpts = int(gridEnhance*np.round(gridspan / lmin))
    #want ngridpts to be an even number so that the noise array can be constructed with same size 
    if np.mod(ngridpts, 2):
        ngridpts += 1

    xvals = np.linspace(0,np.floor(gridspan),ngridpts)
    yvals = xvals



    #Define a falling spectrum 
    Vo = np.logspace(np.log10(p0), np.log10(p1), nfreqs)





    #--------------------------------------------------
    #Final wave values in array of size [N x N x t]. To avoid nested for loops, broadcast all grid and time arrays to this size 

    


    #grid values for each time
    xvals2 = np.repeat(xvals[:,np.newaxis],ngridpts,axis=1)
    yvals2 = xvals2.T

    #This step is a real hog for large arrays
    xvals2 = np.repeat(xvals2[:,:,np.newaxis],ntsteps,axis=2)
    yvals2 = np.repeat(yvals2[:,:,np.newaxis],ntsteps,axis=2)

    tvals2 = np.repeat(tvals[np.newaxis,:],ngridpts, axis=0)
    tvals2 = np.repeat(tvals2[np.newaxis,:,:],ngridpts, axis=0)







    #Calculate the propagating wave
    #print(ntsteps * ngridpts**2)
    V = np.zeros((ngridpts, ngridpts, ntsteps))
    for w in range(nfreqs):
        V += Vo[w]*np.cos(wvals[w]*tvals2 - (kx[w]*xvals2 + ky[w]*yvals2))



    #add in noise
    if not no_noise:
        np.random.seed(0)
        noise = generate_perlin_noise_3d((ngridpts, ngridpts, ntsteps), (int(np.ceil(ngridpts/2)), int(np.ceil(ngridpts/2)), int(np.ceil(ntsteps/2))), tileable=(False, False, False))
        #noise = generate_perlin_noise_3d((ngridpts, ngridpts, ntsteps), (ngridpts/4, ngridpts/4, ntsteps/4), tileable=(False, False, False))
        V += noise




    if animation:
        #Watch the animated wave move! 
        #----Note that the x, y values are indices, not physical distance
        import matplotlib.animation as animation

        V2 = np.transpose(V)  #Need shape [time, xgrid, ygrid]

        fig = plt.figure()
        images = [
            [plt.imshow(layer, cmap='turbo', interpolation='lanczos', animated=True)]
            for layer in V2
        ]
        animation_3d = animation.ArtistAnimation(fig, images, interval=10, blit=True)
        plt.show()





    vals = {'lambda_min':lmin, 'lambda_max':lmax,
            'freqs':fvals,
            'kx':kx,
            'ky':ky,
            'angle':angle,
            'vphase':Vph}

    return(V, xvals, yvals, tvals, vals)



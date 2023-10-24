
"""
Place an E-field antenna system in a wave field of type:
    V(r,t) = sum_i [ Vo_i * cos(w_i*t - k_i*r)]
and measure the probe response (voltage and E-field).
Outputs artificial waveforms that can be used to test various routines
(e.g. interferometry_fvsk_spec_test.py)


NOTE: if wavelength is on order of or shorter than boom length then opposite probes 
won't necessarily be 180 deg out of phase for a monochromatic plane wave

NOTE: phase identification can be difficult if:

1) the wave is propagating nearly along 
one of the interferometry axes (e.g. if along the y-hat' direction then V32 and V14 will
measure very little signal, meaning the waveform can be noise-dominated)

2) if wavelength is on the order or less than boom length than results can't be trusted.


Coordinate system used

                       (y-hat)
                         V3
                      /  |  \             
           (y-hat') x    |     x (x-hat')         
                  /      |       \        
               V2--------O--------V1 (x-hat)    
                  \      |       /
                    x    |     x
                      \  |  /
                         V4
          

x-points represent the centers of potential of the interferometry diagonals          

Coordinate system (system of input test wave)
x-hat --> E12 = V1 - V2 direction (positive to right)
y-hat --> E34 = V3 - V4 direction (positive upwards)

This code outputs the spectrum of k-values in the kx' and ky' directions, where
x'-hat --> Ex' = V1V3x - V4V2x (45 deg inclined from x-hat)
y'-hat --> Ey' = V3V2x - V1V4x (45 deg inclined from y-hat)


"""
import numpy as np 

import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/Endurance/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/signal_analysis/')
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/antenna_analysis/')
import plot_spectrogram as ps
import matplotlib.pyplot as plt
import plane_wave_2d as pw2
from plane_wave_2d import plane_wave_2d as pw2
import correlation_analysis as ca

#Set probe length (probe to sc body). The antenna system will be centered in the wave field
boomL = 5  #length of each Endurance foldout boom



#Determine Vphase based on desired wavelength of plane wave 
l = 50   # desired wavelength meters
f = 20  # freq in Hz
Vph = f*l
gridspan = 2*l
snapshot = 0

#Get 2D waveform that varies with time 
#Grid only needs to be a bit bigger than the tip-tip boom length. 
V,xg,yg,tg,vals = pw2(gridspan=gridspan, #meters
                    f0=f,f1=60,nfreqs=3,
                    angle=165, #propagation angle relative to x-axis
                    Vph=Vph, #phase velocity in m/s
                    snapshot=snapshot,
                    animation=0,
                    timeMultiplier=50,
                    no_noise=0,
                    timeEnhance=16,
                    gridEnhance=16)




#------------------------------------------------------
#Place antennas in the grid and measure the voltage 

#boom length in xg, yg indice units 
grid_per_meter = len(xg) / gridspan
boomLgrid = boomL * grid_per_meter  #boom length as measured in grid spacings (not meters)

posC = [int(len(xg)/2), int(len(yg)/2)]  #center of grid
posT = [posC[0], int(posC[1] + (boomLgrid/2))]     #top probe
posB = [posC[0], int(posC[1] - (boomLgrid/2))]     #bottom probe
posL = [int(posC[0] - (boomLgrid/2)), posC[1]]     #left probe
posR = [int(posC[0] + (boomLgrid/2)), posC[1]]     #right probe


#Extract time-varying voltage for each probe
#V1 = Vp1 - Vsc 
#Apply Hanning window so I can repeat the waveforms to make them larger for a coherence analysis.

Vsc = V[posC[0],posC[1]]
if snapshot != 1:
    hanning = np.hanning(len(Vsc))
    V1 = hanning * (V[posR[0],posR[1],:] - Vsc)
    V2 = hanning * (V[posL[0],posL[1],:] - Vsc)
    V3 = hanning * (V[posT[0],posT[1],:] - Vsc)
    V4 = hanning * (V[posB[0],posB[1],:] - Vsc)
else:
    V1 = V[posR[0],posR[1],:] - Vsc
    V2 = V[posL[0],posL[1],:] - Vsc
    V3 = V[posT[0],posT[1],:] - Vsc
    V4 = V[posB[0],posB[1],:] - Vsc



#Now define the probe pairs that constitute the interferometer
E14 = V1 - V4 
E32 = V3 - V2 

E13 = V1 - V3 
E42 = V4 - V2

E12 = V1 - V2
E34 = V3 - V4

#Plot FFT

#plt.plot(tg,E13)
#plt.xlim(0,0.05)
#ft = np.fft.fft(E13)
#timestep = tg[1]-tg[0]
#freqs = np.fft.fftfreq(ft.shape[-1], d=timestep)
#plt.plot(freqs,np.abs(ft))
#plt.xlim(0,1000)



fig,axs = plt.subplots(3)
axs[0].plot(tg,E12,tg,E34)
axs[1].plot(tg,E14,tg,E32)
axs[2].plot(tg,E13,tg,E42)
for i in range(3): axs[i].set_xlim(2.4,2.8)


"""
#Test snapshot for correct k, lambda
    ax = plt.gca()
    ax.imshow(np.squeeze(V[:,:,0]))
    xticks = np.arange(0,50,20)
    ax.set_xticks(xticks)
    ax.set_yticks(xticks)

    plt.imshow(np.squeeze(V[:,:,0]))
    plt.xticks(xticks)
    plt.yticks(xticks)

    print('here')
"""


#Verify phase velocity 
if not snapshot: 

    fs = 1/(tg[1]-tg[0])

    #NOTE: + sense of phase defined as pointing towards center of potential of first waveform
    ph1,pht1 = ca.phase_cc_timelag_analysis(E13,E42,tg,fs,ccstep=5)
    ph2,pht2 = ca.phase_cc_timelag_analysis(E32,E14,tg,fs,ccstep=5)

    fig,axs = plt.subplots(4)
    axs[0].plot(tg,E13,tg,E42)
    axs[1].plot(pht1,ph1,'.')
    axs[1].plot(pht1,ph1)
    axs[1].set_ylim(np.min(ph1),np.max(ph1))

    axs[2].plot(tg,E32,tg,E14)
    axs[3].plot(pht2,ph2,'.')
    axs[3].plot(pht2,ph2)
    axs[3].set_ylim(np.min(ph2),np.max(ph2))
 
 
#---------------------------------------------
#Once I've identified angles calculate K-vectors
#---------------------------------------------

    phase_xprime = -67 #1342
    phase_yprime = 112  #3214 
    #k-vectors in the primed coord (diagonal interferometer directions)
    kxp = np.radians(phase_xprime)/3.54
    kyp = np.radians(phase_yprime)/3.54
    kmag_measured = np.sqrt(kxp**2 + kyp**2)
    kmag_exp = np.sqrt(vals['kx']**2 + vals['ky']**2)


    #project k-vector in (x',y') coord onto (x,y) 
    angle_rot = np.radians(45) #angle of rotation b/t primed and unprimed coord
                               #i.e. interferometry vs antenna coord
    #(compare to vals['kx'] and vals['ky'])
    kx = kxp*np.cos(angle_rot) + kyp*np.cos(angle_rot + np.radians(90))
    ky = kxp*np.cos(np.radians(90) - angle_rot) + kyp*np.cos(angle_rot)


    #determine measured angle (to compare to angle of test wave). 
    angle_measured = np.degrees(np.arctan2(ky,kx))
    
    wavelength = 2*np.pi/kmag_measured
    wavelength_exp = 2*np.pi/kmag_exp

#Plot hodograms
# --E12 and E34, for long wavelengths, should always be linear and out of phase since they are measuring from the same center of potential
#       They can have other phase for short wavelengths
# --The interferometry pairs will in general measure elliptical polarizations. 
# --In specific case that wave fronts are perp to line b/t center of potentials
# --then they will measure linear
# 

rmax = np.ceil(np.nanmax([E12,E34,E14,E32,E13,E42]))

fig, axs = plt.subplots(3)
axs[0].plot(E12,E34)
#axs[0].set_xlim(-rmax, rmax)
axs[0].set_aspect('equal')
axs[1].plot(E13,E42)
#axs[1].set_xlim(-rmax, rmax)
axs[1].set_aspect('equal')
axs[2].plot(E14,E32)
#axs[2].set_xlim(-rmax, rmax)
axs[2].set_aspect('equal')




import pickle 

wave = {'tvals':tg,
       'wf12':E12,
       'wf34':E34,
       'wf14':E14,
       'wf32':E32,
       'wf13':E13,
       'wf42':E42,
       'vals':vals
       }


goo = pickle.dump(wave, open('/Users/abrenema/Desktop/wavetst.pkl','wb'))



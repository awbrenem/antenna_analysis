
"""
Place an E-field antenna system in a wave field of type:
    V(r,t) = sum_i [ Vo_i * cos(w_i*t - k_i*r)]
and measure the probe response (voltage and E-field).
Outputs artificial waveforms that can be used to test various routines. 

NOTE: if wavelength is on order of or shorter than boom length then opposite probes 
won't necessarily be 180 deg out of phase for a monochromatic plane wave

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


#Set probe length (probe to sc body). The antenna system will be centered in the wave field
boomL = 5  #length of each Endurance foldout boom



#Get 2D waveform that varies with time 
#V,xg,yg,tg,fv,lmin,lmax = pw2(gridspan=30, #meters
#                    f0=2000,f1=4000,nfreqs=100, 
#                    angle=25, #propagation angle relative to x-axis
#                    Vph=1000, #phase velocity in m/s
#                    snapshot=0)


#Get 2D waveform that varies with time 
V,xg,yg,tg,vals = pw2(gridspan=60, #meters
                    f0=100,f1=100,nfreqs=1, 
                    angle=42, #propagation angle relative to x-axis
                    Vph=1000, #phase velocity in m/s
                    snapshot=0,
                    animation=1)




plt.plot(tg, V[0,0,:])


ps.plot_spectrogram(xg, yg, V[:,:,0],vr=[np.min(V),np.max(V)],zscale='linear')


print('h')

#fig, axs = plt.subplots(6)
#ps.plot_spectrogram(xvals, yvals, V[:,:,0],ax=axs[0])
#ps.plot_spectrogram(xvals, yvals, V[:,:,1],ax=axs[1])
#ps.plot_spectrogram(xvals, yvals, V[:,:,2],ax=axs[2])
#ps.plot_spectrogram(xvals, yvals, V[:,:,3],ax=axs[3])
#ps.plot_spectrogram(xvals, yvals, V[:,:,4],ax=axs[4])
#ps.plot_spectrogram(xvals, yvals, V[:,:,5],ax=axs[5])


#ps.plot_spectrogram(xvals, yvals, V[:,:,0])



#------------------------------------------------------
#Place antennas in the grid and measure the voltage 

posC = [int(len(xg)/2), int(len(yg)/2)]  #center of grid
posT = [posC[0], int(posC[1]+boomL)]     #top probe
posB = [posC[0], int(posC[1]-boomL)]     #bottom probe
posL = [int(posC[0]-boomL), posC[1]]     #left probe
posR = [int(posC[0]+boomL), posC[1]]     #right probe


#Extract time-varying voltage for each probe
#V1 = Vp1 - Vsc 
#Apply Hanning window so I can repeat the waveforms to make them larger for a coherence analysis.

Vsc = V[posC[0],posC[1]]
hanning = np.hanning(len(Vsc))
V1 = hanning * (V[posT[0],posT[1],:] - Vsc)
V2 = hanning * (V[posB[0],posB[1],:] - Vsc)
V3 = hanning * (V[posL[0],posL[1],:] - Vsc)
V4 = hanning * (V[posR[0],posR[1],:] - Vsc)





#Now define the probe pairs that constitute the interferometer
E14 = V1 - V4 
E32 = V3 - V2 

E13 = V1 - V3 
E42 = V4 - V2

E12 = V1 - V2
E34 = V3 - V4

plt.plot(tg,V1)
plt.plot(tg,V2)
plt.plot(tg,V3)
plt.plot(tg,V4)

plt.plot(tg,E12)
plt.plot(tg,E34)

plt.plot(tg,E14)
plt.plot(tg,E32)


plt.plot(tg,E13,'.')
plt.plot(tg,E42,'.')
plt.xlim(0.02, 0.022)

#Plot hodograms - should be linear if wavelength >> boom length
plt.plot(V1,V2)
plt.plot(V3,V4)
plt.plot(V1,V3)



#Repeat waveform to build up larger waveform. 
nrepeat = 100


#New times
ln = len(tg)
tmax = np.max(tg)
tg2 = np.zeros(nrepeat * ln)
for r in range(nrepeat):
    tg2[r*ln:(r+1)*ln] = tg + r*tmax


#New E-field components
E14r = np.tile(E14,nrepeat)
E32r = np.tile(E32,nrepeat)
E13r = np.tile(E13,nrepeat)
E42r = np.tile(E42,nrepeat)
E12r = np.tile(E12,nrepeat)
E34r = np.tile(E34,nrepeat)



plt.plot(tg2,E14r)

import pickle 

wave = {'tvals':tg2,
       'wf12':E12r,
       'wf34':E34r,
       'wf14':E14r,
       'wf32':E32r,
       'wf13':E13r,
       'wf42':E42r,
       'vals':vals
       }


goo = pickle.dump(wave, open('/Users/abrenema/Desktop/wavetst.pkl','wb'))



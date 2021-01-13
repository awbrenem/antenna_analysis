;Plots values from rbsp_efw_deconvol_inst_resp.pro
;These characterize the EFW instrument response from signals from the spin plane and spin axis
;booms both with and without AC high-pass filtering.

;**For comparison curves see JBT's document "efw_response", put together pre-launch 
;based on THEMIS EFI. 


;----------QUESTIONS FOR BONNELL 
;1) is the phase vs f correct for the high pass filter? It's different than JBT's
;2) what about the DFB RESPONSE????  This is missing but is included in JBT's plots
;3) Why not calibrating the Vb1 and Vb2 data?



charsz_plot = 0.8             ;character size for plots
charsz_win = 1.2
!p.charsize = charsz_win



rbsp_efw_init
!p.charsize = 1

srate = 16384.

rsheath = 50d6  ;Ohms
rsheathS = 'Rsheath='+string(rsheath/1d6,format='(i3)')+' MOhms'

; Setup kernel for deconvolving EFI response that includes
;   boom response
;   AC high-pass filter (= 1)
;   anti-aliasing Bessel filter response
;   ADC interleaving timing


nelem = 1000.
fmin = 0.001
fmax = 1d6
f = 10d^(indgen(nelem)*alog10(fmax/fmin)/(nelem-1))*fmin
f = [reverse(f),f]


;------------------------------------------
;Determine and plot the boom response 
;Attempting to match Wygant13 Fig10 plot for the 
;curves that represent the input to the IDPU

;Note that this is missing the higher frequency roll-off in the sensor/cable system that occurs as a consequence 
;of the capacitive loading of the preamplifier by the long spin-plane cables. 
;This occurs near 300–400 kHz and defines the upper end of the frequency response 
;for the signals sent to the EMFISIS instrument. -- Highlighted feb 4, 2020
;------------------------------------------

spb_resp = rbsp_efw_boom_response(f, 'SPB', rsheath = rsheath)
axb_resp = rbsp_efw_boom_response(f, 'AXB', rsheath = rsheath)

;Plot boom response (gain and phase)
!p.multi = [0,0,2]
plot,f[0:nelem+1],spb_resp[0:nelem+1],/xlog,xrange=[1d-1,1d6],yrange=[0,1],/ystyle,/xstyle,ytitle='Gain (Vout/Vin)',xtitle='freq(Hz)',title='EFW boom response '+rsheathS,/nodata
oplot,f[0:nelem+1],spb_resp[0:nelem+1],color=250
oplot,f[0:nelem+1],axb_resp[0:nelem+1],color=50
plot,f[0:nelem+1],imaginary(spb_resp[0:nelem+1])/!dtor,/xlog,xrange=[1d-1,1d6],yrange=[-100,5],/ystyle,/xstyle,ytitle='Phase (deg)',xtitle='freq(Hz)',/nodata
oplot,f[0:nelem+1],imaginary(spb_resp[0:nelem+1])/!dtor,color=250
oplot,f[0:nelem+1],imaginary(axb_resp[0:nelem+1])/!dtor,color=50


;------------------------------------------
;Analog filters: Determine and plot the high-pass filter and Bessel anti-aliasing filter
;------------------------------------------

hipass_respDC = 1d
hipass_respAC = rbsp_ac_highpass_response(f)
bessel_resp = rbsp_anti_aliasing_response(f)

!p.multi = [0,2,2]
plot,f[0:nelem+1],hipass_respAC[0:nelem+1],/xlog,xrange=[1d-3,1d4],yrange=[0.1,2],/ylog,/ystyle,/xstyle,ytitle='Gain (Vout/Vin)',xtitle='freq(Hz)',title='hipass AC filter'
plot,f[0:nelem+1],imaginary(hipass_respAC[0:nelem+1])/!dtor,/xlog,xrange=[1d-3,1d4],yrange=[-5,100],/ystyle,/xstyle,ytitle='Phase (deg)',xtitle='freq(Hz)'
plot,f[0:nelem+1],bessel_resp[0:nelem+1],/xlog,xrange=[1d-3,1d4],yrange=[0.1,2],/ylog,/ystyle,/xstyle,ytitle='Gain (Vout/Vin)',xtitle='freq(Hz)',title='5 pole anti-aliasing Bessel filter'
plot,f[0:nelem+1],imaginary(bessel_resp[0:nelem+1])/!dtor,/xlog,xrange=[1d-3,1d4],yrange=[-100,5],/ystyle,/xstyle,ytitle='Phase (deg)',xtitle='freq(Hz)'




;------------------------------------------
;ADC and DFB response
;------------------------------------------

adcresp_12 = rbsp_adc_response(f, 'E12AC') ;;   Calculate the responses of the RBSP DFB ADC chip.
adcresp_34 = rbsp_adc_response(f, 'E34AC')
adcresp_56 = rbsp_adc_response(f, 'E56AC')

!p.multi = [0,2,2]
plot,f[0:nelem+1],adcresp_12[0:nelem+1],/xlog,xrange=[1d-3,1d4],yrange=[0.1,2],/ylog,/ystyle,/xstyle,ytitle='Gain (Vout/Vin)',xtitle='freq(Hz)',title='ADC response E12AC'
plot,f[0:nelem+1],imaginary(adcresp_12[0:nelem+1])/!dtor,/xlog,xrange=[1d-3,1d4],yrange=[-90,10],/ystyle,/xstyle,ytitle='Phase (deg)',xtitle='freq(Hz)'


;****MISSING!!!!!
;***WHAT ABOUT DFB RESPONSE?????



;---------------------------------------------------------
;Overall response (Vout after DFB vs Vin as probe input)
;---------------------------------------------------------

E12_respDC = 1d / (spb_resp * hipass_respDC * bessel_resp * adcresp_12)
E34_respDC = 1d / (spb_resp * hipass_respDC * bessel_resp * adcresp_34)
E56_respDC = 1d / (axb_resp * hipass_respDC * bessel_resp * adcresp_56)


!p.multi = [0,2,2]
plot,f[0:nelem+1],E12_respDC[0:nelem+1],/xlog,xrange=[1d-3,1d4],yrange=[0.1,2],/ylog,/ystyle,/xstyle,ytitle='Gain (Vout/Vin)',xtitle='freq(Hz)',title='E12DC response'
plot,f[0:nelem+1],imaginary(E12_respDC[0:nelem+1])/!dtor,/xlog,xrange=[1d-3,1d4],yrange=[-200,200],/ystyle,/xstyle,ytitle='Phase (deg)',xtitle='freq(Hz)'


;Eb2 (AC-coupled)
;   AC high-pass filter
hipass_respAC = 1d
E12_respAC = 1d / (spb_resp * hipass_respAC * bessel_resp * adcresp_12)
E34_respAC = 1d / (spb_resp * hipass_respAC * bessel_resp * adcresp_34)
E56_respAC = 1d / (axb_resp * hipass_respAC * bessel_resp * adcresp_56)



ind = where(finite(E12_respDC, /nan), nind) & if nind gt 0 then E12_respDC[ind] = 0
ind = where(finite(E34_respDC, /nan), nind) & if nind gt 0 then E34_respDC[ind] = 0
ind = where(finite(E56_respDC, /nan), nind) & if nind gt 0 then E56_respDC[ind] = 0
ind = where(finite(E12_respAC, /nan), nind) & if nind gt 0 then E12_respAC[ind] = 0
ind = where(finite(E34_respAC, /nan), nind) & if nind gt 0 then E34_respAC[ind] = 0
ind = where(finite(E56_respAC, /nan), nind) & if nind gt 0 then E56_respAC[ind] = 0


; Transfer kernel into time domain: take inverse FFT and center
E12_respDC2 = shift((fft(E12_respDC,1)), kernel_length/2) / kernel_length
E34_respDC2 = shift((fft(E34_respDC,1)), kernel_length/2) / kernel_length
E56_respDC2 = shift((fft(E56_respDC,1)), kernel_length/2) / kernel_length

E12_respAC2 = shift((fft(E12_respAC,1)), kernel_length/2) / kernel_length
E34_respAC2 = shift((fft(E34_respAC,1)), kernel_length/2) / kernel_length
E56_respAC2 = shift((fft(E56_respAC,1)), kernel_length/2) / kernel_length


;Mimic Bonnell08 EFI plot in Fig7

!p.multi = [0,0,2]
plot,f[0:511],e12_respDC[0:511],xrange=[10,10000],/xlog,title='EB1 Transfer function (output/input)!CBlack=SPB,Blue=AXB'
oplot,f[0:511],real_part(e12_respDC[0:511]),color=50
;plot,f[0:511],bessel_resp,title='RE bessel resp (soft rolloff filter)',xrange=[10,10000],/xlog
;plot,f[0:511],spb_resp,title='spb resp (antenna coupling to plasma: RC circuit)',xrange=[10,10000],/xlog
;plot,f[0:511],adcresp_12,title='ADC resp',xrange=[10,10000],/xlog,xtitle='freq (Hz)'


!p.multi = [0,0,2]
plot,f[0:511],real_part(bessel_resp),title='RE bessel resp (soft rolloff filter)',xrange=[10,10000],/xlog
plot,f[0:511],imaginary(bessel_resp),title='IM bessel resp (soft rolloff filter)',xrange=[10,10000],/xlog,xtitle='freq (Hz)'





;compare integration test curve (ADC resp only) to the full transfer function.
adc_test = 1/adcresp_12
Vout_Vin_gain = abs(e12_resp/adc_test)
Vout_Vin_phase = imaginary(e12_resp/adc_test)
Vout_Vin_gain = abs(e12_resp)
Vout_Vin_phase = imaginary(e12_resp)
;ratio2 = abs(e12_resp/adc_test)

!p.multi = [0,0,3]
plot,f[0:511],1/(e12_resp[0:511]/adc_test),yrange=[0,1],xrange=[10,10000],/xlog,title='curve1/curve2. Amt values will be off if only ADCgain is accounted for'
plot,f[0:511],Vout_Vin_gain[0:511],yrange=[0,1],xrange=[10,10000],/xlog,title='curve1/curve2. Amt values will be off if only ADCgain is accounted for'
plot,f[0:511],1/Vout_Vin_phase[0:511],yrange=[-100,100],xrange=[10,10000],/xlog,title='curve1/curve2. Amt values will be off if only ADCgain is accounted for'

plot,f[0:511],e12_resp[0:511],xrange=[10,10000],/xlog,title='Transfer function'
plot,f[0:511],adc_test[0:511],xrange=[10,10000],/xlog,title='1/ADC_resp12'
plot,f[0:511],adcresp_12[0:511],xrange=[10,10000],/xlog,title='ADC_resp12'




;;Plot without considering poles in complex plane
!p.multi=[0,0,1]
vals = real_part(e12_resp[0:511])
vals2 = 1/vals
plot,f[0:511],vals,xrange=[10,10000],/xlog
oplot,f[0:511],vals2,color=250;xrange=[10,10000],/xlog

;Take into account poles in complex plane
vals = e12_resp[0:511]
vals2 = 1/vals
plot,f[0:511],vals,xrange=[10,10000],/xlog,title='Transfer function (black)'
oplot,f[0:511],vals2,color=250;xrange=[10,10000],/xlog





!p.multi = [0,0,5]
plot,f[0:511],e12_respDC[0:511],xrange=[10,10000],/xlog,title='EB2 Transfer function'
plot,f[0:511],bessel_resp,title='bessel resp (soft rolloff filter)',xrange=[10,10000],/xlog
plot,f[0:511],spb_resp,title='spb resp (antenna coupling to plasma: RC circuit)',xrange=[10,10000],/xlog
plot,f[0:511],adcresp_12,title='ADC resp',xrange=[10,10000],/xlog
plot,f[0:511],hipass_resp[0:511],xrange=[10,10000],/xlog,title='High pass filter response'



;oplot,e12_resp,color=250
oplot,bessel_resp,color=250

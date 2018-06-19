;Plots the STEREO TDS gain and phase curves as a function of freq 
;from patdsgain.pro

;These gain/phase curves don't take into account the antenna impedence (i.e.
;the coupling to the plasma). It's just the gain of TDS and the preamp.


pro stereo_gain_phase_plot

	tempf = 30.5176
	npoints = 4096
;    tempf=1.0/(deltat*npoints)
    tds_corr=complexarr(npoints)
    mvc_temp=fltarr(npoints)
    tds_corr[0]=complex(1.0,0.0)
    tds_corr[npoints/2]=tds_corr[0]
    ipa=1
    ifilt = 0.  ;value from swtds.pro
	freqs = fltarr(npoints)
	phase = fltarr(npoints)

	gain = fltarr(npoints)
    for ccount=1,npoints/2-1 do begin  $
      f=tempf*ccount   & $
      patdsgain,f,ipa,ifilt,cgaintds,gaintds,phasetds  & $
	  freqs[ccount] = f & $
      tds_corr[ccount]=cgaintds & $
      tds_corr[npoints-ccount]=conj(cgaintds)  & $
	  gain[ccount] = gaintds      & $
	  phase[ccount] = phasetds  & $
    endfor


!p.multi = [0,0,2]
plot,freqs,gain,/xlog,xrange=[10.,1d5],title='STEREO TDS GAIN VS FREQ (kHz) (from stereo_gain_phase_plot.pro)'
plot,freqs,phase,/xlog,xrange=[10.,1d5],title='STEREO TDS PHASE VS FREQ (kHz) (from stereo_gain_phase_plot.pro)'




end
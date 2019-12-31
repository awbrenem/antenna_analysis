;Creates comparison plots of the single-ended probe potentials for
;RBSP EFW. This includes a cross-correlation analysis of the phase
;difference for each antenna combination

;Also determines times of wake effects (via Sheng's code, which he says
;works very well), and compares ratio of Eu to Ev. Flags times that are not
;too close to unity.

;Note on phase differences b/t antenna combinations
;(1) IN THE PRESENSE OF EXTERNAL E-FIELD: Opposite antenna pairs should be 180 deg different (V1-V2 and V3-V4), while
;   (V1-V3) and (V2-V4) and (V1-V4) should be 90 or 270.
;(2) IN THE PRESENSE OF DENSITY FLUCTUATIONS: Antennas can all have the same phase.
;(3) EXTERNAL E-FIELD + DENSITY --> various phase differences depending on relative
;   strength of the two.

;Note on warped sinusoids: This can mean the presense of a wake effect. Electrons
;and ions can have different thermal temps, allowing one to fill in wake faster than
;the other. This results in an artificial Efield.


;Keywords:
;   cc_stepsize --> number of spinperiods to calculate the cross-correlations over.
;              Defaults to 10.
;   detrend_time --> time (sec) to detrend over

;; 0.020 Hz  -> 50 sec period


pro efw_probe_potentials_test,probe=probe,noplot=noplot,detrend_time=dettime,$
  cc_stepsize=ccstep

  tplot_options,'xmargin',[20.,16.] & tplot_options,'ymargin',[3,4]
  tplot_options,'xticklen',0.08 & tplot_options,'yticklen',0.02
  tplot_options,'xthick',2 & tplot_options,'ythick',2 & tplot_options,'labflag',-1


  pairs = ['12','13','14','23','24','34']
  npairs = n_elements(pairs)-1


  ;-------------------------------------------------------------------------
  ;Keeping this number low (e.g. 1 or 2) you can still see if the phases are in sync
  ;during times when the probes are fluctuating by density fluctuations.
  ;These values tend to be fairly choppy. Later on I smooth these values to reduce this
  ;choppiness.
  ;Can also set to a higher number (e.g. 10 or 20). The relative phases will tend to
  ;go to zero when there are density fluctuations. No need to detrend these.
  ;  if ~KEYWORD_SET(ccstep) then ccstep = 10.  ;# of spinperiods to calculate cc over

  if ~KEYWORD_SET(ccstep) then ccstep = 4.  ;# of spinperiods to calculate cc over
  ccstep_detrend_spinperiods = 20.
  ;-------------------------------------------------------------------------

  ;Number of spin periods to calculate the Efield envelope over
  n_spinperiods_envelope = 5.



  tr = time_string(timerange()) & tr = tr[0]
  date = strmid(tr,0,10)
  rbx = 'rbsp'+probe+'_'


  ;; get spinperiod
  rbsp_efw_hsk_crib,date,probe
  get_data,rbx+'efw_hsk_idpu_eng_ACS_PERIOD',data=dat
  spinperiod = median(dat.y)


  if ~KEYWORD_SET(dettime) then dettime = 6.*spinperiod



  ;---------------------------------------------
  ;Load L3 data to get spinfit values and flags
  ;---------------------------------------------

  rbsp_load_efw_waveform_l3,probe=probe
  split_vec,rbx+'efw_flags_charging_bias_eclipse',suffix=['_charging','_extremecharging','_bias','_eclipse']

  options,rbx+'efw_flags_charging_bias_eclipse_charging','ytitle',rbx+'!Ccharging!Cflag'
  options,rbx+'efw_flags_charging_bias_eclipse_bias','ytitle',rbx+'!Cbias!Cflag'
  options,rbx+'efw_flags_charging_bias_eclipse_eclipse','ytitle',rbx+'!Ceclipse!Cflag'
  options,rbx+'efw_flags_charging_bias_eclipse_extremecharging','ytitle',rbx+'!Cextremecharging!Cflag'
  ylim,[rbx+'efw_flags_all',rbx+'efw_flags_charging_bias_eclipse',$
  rbx+'efw_flags_charging_bias_eclipse_*'],0,2
  ylim,rbx+'efw_density',8,10000,1
  ylim,rbx+'efw_efield_inertial_frame_mgse',-10,10



  ;--------------------------------------------------
  ;CHECK SINGLE-ENDED MEASUREMENTS
  ;--------------------------------------------------

  rbsp_load_efw_waveform,probe=probe,type='calibrated',datatype='vsvy'
  split_vec,rbx+'efw_vsvy',suffix='_'+['v1','v2','v3','v4','v5','v6']



  ;; get samplerate
  get_data,rbx+'efw_vsvy_v1',data=dat
  sr = rbsp_sample_rate(dat.x,out_med_avg=medavg)
  sr = medavg[0]

  ;-----------------------------------------------------------------------
  ;Create bandpassed versions of the single-ended potentials. These will be
  ;used to calculate the cross-correlation. This is much cleaner than using
  ;detrend.

  ;smooth high freq fluctuations
  rbsp_detrend,rbx+'efw_vsvy_v?',2
  ;get rid of DC offset
  rbsp_detrend,rbx+'efw_vsvy_v?_smoothed',dettime

  copy_data,rbx+'efw_vsvy_v1_smoothed_detrend',rbx+'v1_bp'
  copy_data,rbx+'efw_vsvy_v2_smoothed_detrend',rbx+'v2_bp'
  copy_data,rbx+'efw_vsvy_v3_smoothed_detrend',rbx+'v3_bp'
  copy_data,rbx+'efw_vsvy_v4_smoothed_detrend',rbx+'v4_bp'




  for i=0,npairs do store_data,rbx+'comb'+pairs[i]+'_bp',data=[rbx+'v'+strmid(pairs[i],0,1)+'_bp',rbx+'v'+strmid(pairs[i],1,1)+'_bp']
  for i=0,npairs do options,rbx+'comb'+pairs[i]+'_bp','ytitle','RBSP'+probe+'!C'+'V'+strmid(pairs[i],0,1)+',V'+strmid(pairs[i],1,1)+'!Cdetrend'
  options,[rbx+'comb??_bp'],'colors',[0,250]


  for i=0,npairs do store_data,rbx+'comb'+pairs[i],data=[rbx+'efw_vsvy_v'+strmid(pairs[i],0,1),rbx+'efw_vsvy_v'+strmid(pairs[i],1,1)]
  for i=0,npairs do store_data,rbx+'comb'+pairs[i]+'d',data=[rbx+'efw_vsvy_v'+strmid(pairs[i],0,1),rbx+'efw_vsvy_v'+strmid(pairs[i],1,1)]+'_smoothed_detrend'


  for q=0,npairs do begin

    vat = strmid(pairs[q],0,1)
    vbt = strmid(pairs[q],1,1)

    ;density proxy
    get_data,rbx+'efw_vsvy_v'+vat,data=v1
    get_data,rbx+'efw_vsvy_v'+vbt,data=v2
    store_data,rbx+'density'+pairs[q]+'_proxy',data={x:v1.x,y:(v1.y+v2.y)/2.}

    rbsp_efw_density_fit_from_uh_line,rbx+'density'+pairs[q]+'_proxy',probe,$
                                      newname=rbx+'density'+pairs[q]+'_actual',$
                                      dmin=0.1,$
                                      dmax=20000.,$
                                      setval=-1.e31

    ylim,rbx+'density'+pairs[q]+'_actual',1,50000,1


    ;    rbsp_detrend,rbx+'density'+pairs[q]+'_proxy',2
    ;    rbsp_detrend,rbx+'density'+pairs[q]+'_proxy_smoothed',60.
    ;    options,rbx+'density'+pairs[q]+'_proxy','ytitle','"Density'+pairs[q]+'proxy"!Cdetrended'
    ;    options,rbx+'density'+pairs[q]+'_proxy_smoothed_detrend','ytitle','"Density'+pairs[q]+'proxy"!Cdetrended'
    ;    ylim,rbx+'density'+pairs[q]+'_proxy_smoothed_detrend',-2,2

  endfor




  ;------------------------------------------------------------
  ;Cross-correlation to determine lag b/t various antennas
  ;------------------------------------------------------------

  get_data,rbx+'v1_bp',data=v1
  get_data,rbx+'v2_bp',data=v2
  get_data,rbx+'v3_bp',data=v3
  get_data,rbx+'v4_bp',data=v4

  lagmax = spinperiod * sr      ;number of data points in a single spinperiod
  lag = indgen(200)*lagmax/(199.)

  stepsize = spinperiod*ccstep
  nsteppts = stepsize * sr

  nsteps = floor((max(v1.x,/nan) - min(v1.x,/nan))/stepsize)


  bb = 0L
  uu = nsteppts
  phase12 = fltarr(nsteps) & phase13 = fltarr(nsteps) & phase14 = fltarr(nsteps)
  phase23 = fltarr(nsteps) & phase24 = fltarr(nsteps) & phase34 = fltarr(nsteps)
  phaseT = fltarr(nsteps)

  for qq=0L,nsteps-1 do begin

    v1t = v1.y[bb:uu] & v2t = v2.y[bb:uu] & v3t = v3.y[bb:uu] & v4t = v4.y[bb:uu]

    v1tx = v1.x[bb:uu]

    result12 = c_correlate(v1t, v2t, lag)
    result13 = c_correlate(v1t, v3t, lag)
    result14 = c_correlate(v1t, v4t, lag)
    result23 = c_correlate(v2t, v3t, lag)
    result24 = c_correlate(v2t, v4t, lag)
    result34 = c_correlate(v3t, v4t, lag)

    maxv12 = max(result12,wh12) & maxv13 = max(result13,wh13) & maxv14 = max(result14,wh14)
    maxv23 = max(result23,wh23) & maxv24 = max(result24,wh24) & maxv34 = max(result34,wh34)

    ;Peak delay from cross-correlation
    delay12T = lag[wh12]/sr
    delay13T = lag[wh13]/sr
    delay14T = lag[wh14]/sr
    delay23T = lag[wh23]/sr
    delay24T = lag[wh24]/sr
    delay34T = lag[wh34]/sr

    ;; plot,lag/sr,result,ytitle='Cross-correlation',xtitle='lag (sec)'
    ;; oplot,[lag[wh]/sr,lag[wh]/sr],[result[wh],result[wh]],psym=4,symsize=3

    ;delay time as fraction of spinperiod
    frac_spinperiod12 = 100*delay12T/spinperiod
    frac_spinperiod13 = 100*delay13T/spinperiod
    frac_spinperiod14 = 100*delay14T/spinperiod
    frac_spinperiod23 = 100*delay23T/spinperiod
    frac_spinperiod24 = 100*delay24T/spinperiod
    frac_spinperiod34 = 100*delay34T/spinperiod

    ;turn the delay time into a phase
    phase12[qq] = 360.*frac_spinperiod12/100.
    phase13[qq] = 360.*frac_spinperiod13/100.
    phase14[qq] = 360.*frac_spinperiod14/100.
    phase23[qq] = 360.*frac_spinperiod23/100.
    phase24[qq] = 360.*frac_spinperiod24/100.
    phase34[qq] = 360.*frac_spinperiod34/100.

    phaseT[qq] = (v1tx[n_elements(v1tx)-1] + v1tx[0])/2.

    bb = bb + nsteppts
    uu = uu + nsteppts
  endfor


  store_data,rbx+'phase_v1v2',data={x:phaseT,y:phase12}
  store_data,rbx+'phase_v1v3',data={x:phaseT,y:phase13}
  store_data,rbx+'phase_v1v4',data={x:phaseT,y:phase14}
  store_data,rbx+'phase_v2v3',data={x:phaseT,y:phase23}
  store_data,rbx+'phase_v2v4',data={x:phaseT,y:phase24}
  store_data,rbx+'phase_v3v4',data={x:phaseT,y:phase34}

  ;Smooth out phase values.
  rbsp_detrend,rbx+'phase_v?v?',10.8*ccstep_detrend_spinperiods

  store_data,'line0',data={x:phaseT,y:replicate(0.,n_elements(phaseT))}
  store_data,'line90',data={x:phaseT,y:replicate(90.,n_elements(phaseT))}
  store_data,'line180',data={x:phaseT,y:replicate(180.,n_elements(phaseT))}
  store_data,'line270',data={x:phaseT,y:replicate(270.,n_elements(phaseT))}
  store_data,'line360',data={x:phaseT,y:replicate(360.,n_elements(phaseT))}
  store_data,'line0red',data={x:phaseT,y:replicate(0.,n_elements(phaseT))}
  store_data,'line90red',data={x:phaseT,y:replicate(90.,n_elements(phaseT))}
  store_data,'line180red',data={x:phaseT,y:replicate(180.,n_elements(phaseT))}
  store_data,'line270red',data={x:phaseT,y:replicate(270.,n_elements(phaseT))}
  store_data,'line360red',data={x:phaseT,y:replicate(360.,n_elements(phaseT))}
  options,['line90red','line180red','line270red','line360red'],'colors',250
  options,['line0','line90','line180','line270','line360'],'colors',50
  options,['line90red','line180red','line270red','line360red'],'thick',4
  options,['line0','line90','line180','line270','line360'],'thick',1

  store_data,rbx+'phase_v1v2c',data=['line0','line90','line180red','line270','line360',rbx+'phase_v1v2_smoothed']
  store_data,rbx+'phase_v1v3c',data=['line0','line90','line180','line270red','line360',rbx+'phase_v1v3_smoothed']
  store_data,rbx+'phase_v1v4c',data=['line0','line90red','line180','line270red','line360',rbx+'phase_v1v4_smoothed']
  store_data,rbx+'phase_v2v3c',data=['line0','line90red','line180','line270','line360',rbx+'phase_v2v3_smoothed']
  store_data,rbx+'phase_v2v4c',data=['line0','line90','line180','line270red','line360',rbx+'phase_v2v4_smoothed']
  store_data,rbx+'phase_v3v4c',data=['line0','line90','line180red','line270','line360',rbx+'phase_v3v4_smoothed']
  ylim,rbx+['phase_v?v?c'],-10,380
  options,rbx+'phase_v?v?','thick',2

  path = '~/Desktop/code/Aaron/RBSP/EFW_data_check/'
  fn = 'antenna_phase_test_RBSP'+probe+'_'+strmid(date,0,4)+strmid(date,5,2)+strmid(date,8,2)


  tplot_options,'title',''
  !p.charsize = 1.
  tplot,[rbx+'efw_efield_inertial_frame_mgse',$
    rbx+'density12_actual',$
    rbx+'efw_flags_charging_bias_eclipse_charging',$
    rbx+'efw_flags_charging_bias_eclipse_bias',$
    rbx+'efw_flags_charging_bias_eclipse_eclipse',$
    rbx+'comb12_bp',$
    rbx+'comb34_bp',$
    rbx+'phase_v1v2c',$
    rbx+'phase_v3v4c']


stop


  ;-------------------------------------------
  ;Check to see if we have wake effects (Sheng's code)
  ;-------------------------------------------

  get_data,rbx+'efw_vsvy_v1',data=v1
  get_data,rbx+'efw_vsvy_v2',data=v2
  get_data,rbx+'efw_vsvy_v3',data=v3
  get_data,rbx+'efw_vsvy_v4',data=v4
  store_data,rbx+'e12_rough',data={x:v1.x,y:10.*(v1.y-v2.y)}
  store_data,rbx+'e13_rough',data={x:v1.x,y:10.*(v1.y-v3.y)}
  store_data,rbx+'e14_rough',data={x:v1.x,y:10.*(v1.y-v4.y)}
  store_data,rbx+'e23_rough',data={x:v1.x,y:10.*(v2.y-v3.y)}
  store_data,rbx+'e24_rough',data={x:v1.x,y:10.*(v2.y-v4.y)}
  store_data,rbx+'e34_rough',data={x:v1.x,y:10.*(v3.y-v4.y)}



  sr = 1/(v1.x[1]-v1.x[0])
;  tplot_options,'title','wave effects test'

  for i=0,npairs do stplot_analysis_spinplane_efield,rbx+'e'+pairs[i]+'_rough',spinrate=sr
  for i=0,npairs do options,rbx+'e'+pairs[i]+'_rough0_flag','ytitle','e12!Cwake flag'
  for i=0,npairs do options,rbx+'e'+pairs[i]+'_rough0_flag','ytitle','e12!Cwake flag'
;  tplot,['rbspa_efw_vsvy_v1','e120','e120_mat','e120_amp','e120_flag']



  date_flags = strmid(tr,0,4)+strmid(tr,5,2)+strmid(tr,8,2)




;--------------------------------------------------------------------
;  Another useful check is to compare the amplitude mismatch between the sine
;  waves of Eu and Ev, which detects more bad data in addition to wakes.
;--------------------------------------------------------------------


rbsp_detrend,rbx+['e12_rough','e34_rough'],spinperiod*4.
store_data,'euv_comb',data=rbx+['e12_rough_detrend','e34_rough_detrend']
options,rbx+'euv_rough_comb','colors',[0,250]

get_data,rbx+'e12_rough_detrend',data=e12det
get_data,rbx+'e34_rough_detrend',data=e34det
sr = 1/(e12det.x[1] - e12det.x[0])   ;S/sec
win_width = n_spinperiods_envelope*spinperiod*sr  ;#points containing a single sine wave
t0 = min(e12det.x)
t1 = max(e12det.x)

envelope,(e12det.y),win_width,mn12,mx12,ind12
ttmp = (t1-t0)*indgen(n_elements(ind12))/(n_elements(ind12)-1) + t0

store_data,rbx+'e12_rough_envelopemx',ttmp,mx12
store_data,rbx+'e12_rough_envelopemn',ttmp,mn12
options,rbx+'e12_rough_envelope??','thick',2

envelope,(e34det.y),win_width,mn34,mx34,ind34
store_data,rbx+'e34_rough_envelopemx',ttmp,mx34
store_data,rbx+'e34_rough_envelopemn',ttmp,mn34
options,rbx+'e34_rough_envelope??','thick',2
options,rbx+'e34_rough_envelope??','color',250
envelope_diffmx = abs(mx12 - mx34)
envelope_diffmn = abs(mn12 - mn34)
envelope_divmx = abs(mx12/mx34) > abs(mx34/mx12)  ;catches both large e12/e34 and large e34/e12
envelope_divmn = abs(mn12/mn34) > abs(mn34/mn12)


;Create the flag value when the envelopes of e12-e34 exceed a certain value
envelope_diff_max = 5. ;mV/m
envelope_div_max = 2. ;ratio
envelope_diff_flag = bytarr(n_elements(mx12))
envelope_div_flag = bytarr(n_elements(mx12))
goo = where((envelope_diffmx ge envelope_diff_max) or (envelope_diffmn ge envelope_diff_max))
if goo[0] ne -1 then envelope_diff_flag[goo] = 1
goo = where(envelope_diffmn ge envelope_diff_max)
if goo[0] ne -1 then envelope_diff_flag[goo] = 1

;Find when ratio (e12/e34) is too high
goo = where((envelope_divmx ge envelope_div_max) or (envelope_divmn ge envelope_div_max))
if goo[0] ne -1 then envelope_div_flag[goo] = 1
goo = where(envelope_divmn ge envelope_div_max)
if goo[0] ne -1 then envelope_div_flag[goo] = 1
;Find when ratio (e12/e34) is too low (e.g. e34 >> e12)
goo = where((envelope_divmx le 1/envelope_div_max) or (envelope_divmn le 1/envelope_div_max))
if goo[0] ne -1 then envelope_div_flag[goo] = 1
goo = where(envelope_divmn ge envelope_div_max)
if goo[0] ne -1 then envelope_div_flag[goo] = 1


store_data,rbx+'envelope_diff_flag',ttmp,envelope_diff_flag
store_data,rbx+'envelope_div_flag',ttmp,envelope_div_flag
ylim,rbx+['envelope_diff_flag','envelope_div_flag'],0,2
options,rbx+'envelope_diff_flag','ytitle','e12-e34!Cenvelope_diff_flag'
options,rbx+'envelope_div_flag','ytitle','e12/e34!Cenvelope!Cflag'


store_data,rbx+'e12_rough_envelopemxmn_test',data=rbx+['e12_rough_envelopemx','e12_rough_envelopemn','e12_detrend']
store_data,rbx+'e34_rough_envelopemxmn_test',data=rbx+['e34_rough_envelopemx','e34_rough_envelopemn','e34_detrend']
;store_data,'e12_envelopemn_test',data=['e12_envelopemn','e12_detrend']
;store_data,'e34_envelopemn_test',data=['e34_envelopemn','e34_detrend']
options,rbx+'e??_rough_envelope??_test','colors',[250,0]
store_data,rbx+'e12_e34_envelopemx_comp',data=rbx+['e12_rough_envelopemx','e34_rough_envelopemx']
store_data,rbx+'e12_e34_envelopemn_comp',data=rbx+['e12_rough_envelopemn','e34_rough_envelopemn']
store_data,rbx+'e12_e34_envelopemx_diff',ttmp,envelope_diffmx
store_data,rbx+'e12_e34_envelopemn_diff',ttmp,envelope_diffmn
store_data,rbx+'e12_e34_envelopemx_div',ttmp,envelope_divmx
store_data,rbx+'e12_e34_envelopemn_div',ttmp,envelope_divmn
options,rbx+'e12_e34_envelope??_comp','colors',[0,250]

store_data,'diffline',ttmp,replicate(envelope_diff_max,n_elements(ttmp))
store_data,'divline',ttmp,replicate(envelope_div_max,n_elements(ttmp))
options,['diffline','divline'],'linestyle',2
store_data,rbx+'envelopemx_diff_comb',data=[rbx+'e12_e34_envelopemx_diff',rbx+'e12_e34_envelopemn_diff','diffline']
store_data,rbx+'envelopemx_div_comb',data=[rbx+'e12_e34_envelopemx_div',rbx+'e12_e34_envelopemn_div','divline']


ylim,[rbx+'e12_e34_envelopemx_div',rbx+'e12_e34_envelopemn_div',rbx+'envelopemx_div_comb'],0.8,10,1

;Plot envelopes for Eu/Ev ratio determination
tplot,[rbx+'efw_efield_inertial_frame_mgse',$
  rbx+'density12_actual',$
  rbx+'efw_flags_charging_bias_eclipse_charging',$
  rbx+'efw_flags_charging_bias_eclipse_bias',$
  rbx+'efw_flags_charging_bias_eclipse_eclipse',$
  rbx+'envelope_diff_flag',rbx+'envelope_div_flag',rbx+'e120_flag',rbx+'e340_flag',$
  rbx+'e12_envelopemxmn_test',rbx+'e34_envelopemxmn_test',$
  rbx+'e12_e34_envelopemx_comp',rbx+'e12_e34_envelopemn_comp',$
  rbx+'envelopemx_diff_comb',rbx+'envelopemx_div_comb']
;tplot,['envelope_diff_flag','e12_e34_envelopemx_comp','e12_e34_envelopemn_comp','e12_e34_envelopemx_diff','e12_e34_envelopemn_diff']
;tplot,['envelope_div_flag','e12_e34_envelopemx_comp','e12_e34_envelopemn_comp','e12_e34_envelopemx_div','e12_e34_envelopemn_div']
;tplot,['envelope_diff_flag','envelope_div_flag','e12_envelopemn_test','e34_envelopemn_test','e12_e34_envelopemn_comp','e12_e34_envelopemn_diff']


;Create ASCII files of these flag values
;  fn_envelopeflag = '~/Desktop/e12_e34_envelope_ratio_flag_RBSP'+probe+'_'+date_flags+'.txt'
;  get_data,'envelope_diff_flag',data=ef
;  get_data,'envelope_div_flag',data=ev
;  openw,lun1,fn_envelopeflag,/get_lun
;  for i=0,n_elements(ef.x)-1 do printf,lun1,time_string(ef.x[i]),ef.y[i],ev.y[i]
;  close,lun1
;  free_lun,lun1


;stop
;--------------------------------------------------------------------
;FINAL PLOTS
;--------------------------------------------------------------------

    ;Smoothed and detrended versions
    tplot,[rbx+'efw_efield_inertial_frame_mgse',$
      rbx+'density12_actual',$
      rbx+'density34_actual',$
      rbx+'efw_flags_all',$
      rbx+'efw_flags_charging_bias_eclipse_charging',$
      rbx+'efw_flags_charging_bias_eclipse_bias',$
      rbx+'efw_flags_charging_bias_eclipse_eclipse',$
      rbx+'envelope_div_flag',$
      rbx+'comb12_bp',rbx+'comb34_bp',$
      rbx+'e120_flag',rbx+'e340_flag',$
      rbx+'phase_v1v2c',rbx+'phase_v3v4c']
    stop
    ;--------------------------------
    ;Non smoothed or detrended versions

    tplot,[rbx+'efw_efield_inertial_frame_mgse',$
        rbx+'efw_vsvy_v1',$
        rbx+'efw_vsvy_v2',$
        rbx+'efw_vsvy_v3',$
        rbx+'efw_vsvy_v4',$
        rbx+'density12_actual',$
        rbx+'density34_actual',$
        rbx+'envelope_div_flag',$
        rbx+'e12_rough0_flag',rbx+'e34_rough0_flag',$
        rbx+'efw_flags_all',$
        rbx+'efw_flags_charging_bias_eclipse_charging',$
        rbx+'efw_flags_charging_bias_eclipse_bias',$
        rbx+'efw_flags_charging_bias_eclipse_eclipse',$
        rbx+'phase_v1v2c',rbx+'phase_v3v4c']

    stop
    ;--------------------------------

  store_data,tnames(['*hsk*','*rh*','*yh*','*yl*','*rl*','*BEB*','*DFB*']),/del;,$
  store_data,['rbsp?_efw_density','rbsp?_efw_Vavg','rbsp?_efw_','*smoothed_smoothed*','*proxy*','rbsp?_efw_vsvy'],/del


  tplot_save,'*',filename='~/Desktop/RBSP'+probe+'_'+date_flags+'_probe_potentials_test'

stop

end

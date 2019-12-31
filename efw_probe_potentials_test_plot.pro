;Plot output .tplot files from efw_probe_potentials_test.pro.
;This is mainly used so that people interested in EFW data quality
;can plot this output.

;Neethal's wondering whether the Efield fluctuations (e.g. near 11UT) are caused by the
;density fluctuations near PP crossings.

;-RBSPb and RBSPa very close together
;-Just post noon at apogee.
;-Quiet day with no charging.

timespan,'2017-08-16/08:00',6,/hours


;Restore .tplot file with the variables
probe = 'b'
tplot_restore,filenames='/Users/aaronbreneman/Desktop/Research/OTHER/Stuff_for_other_people/Thomas_Neethal/RBSPb_20170816_probe_potentials_test.tplot'

;Make plots prettier
rbsp_efw_init
rbspx = 'rbsp'+probe+'_'


rbsp_detrend,'rbsp'+probe+'_efw_vsvy_v?',60.*10.
ylim,'rbsp'+probe+'_efw_vsvy_v?'+'_detrend',-0.4,0.4




;Compare single-ended measurements
tplot,'rbsp'+probe+'_efw_vsvy_v?'+'_detrend'
;NOTE: V3 weaker than others


;Compare all antenna pair combos
tplot,rbspx+'comb??'
;NOTE: seems that V1,V2 is the pair with best comparable signals


;Compare electric fields
tplot,rbspx+'e??_rough'
;NOTE: E34 noisy as expected

;Compare density measurements
tplot,rbspx+'density??_actual'


;Compare E12 with density12
ylim,rbspx+'efw_efield_inertial_frame_mgse',-1,2
tplot,rbspx+['efw_efield_inertial_frame_mgse','e12_rough','density12_actual']



tplot,rbspx+['envelope_diff_flag',$
  'envelope_div_flag',$
  'e12_rough0_flag','e34_rough0_flag',$
  'e12_rough_envelopemxmn_test',$
  'e34_rough_envelopemxmn_test',$
  'e12_e34_envelopemx_comp',$
  'e12_e34_envelopemn_comp',$
  'envelopemx_diff_comb',$
  'envelopemx_div_comb']



tplot,rbspx+['envelope_div_flag',$ ;flag for when E12/E34 ratio varies too far from unity
  'e12_rough0_flag',$                ;flag indicating probable wake effects
  'comb12_bp',$  ;smoothed and detrended V1 and V2
  'comb34_bp',$  ;smoothed and detrended V3 and V4
  'comb13_bp',$  ;smoothed and detrended V1 and V3
  'density12_proxy_smoothed_detrend',$ ;density proxy (V1+V2)/2
  'phase_v1v2c',$ ;phase difference of V1 and V2
  'phase_v3v4c',$ ;phase difference of V3 and V4
  'phase_v1v3c',$ ;phase difference of V1 and V3
  'phase_v2v3c',$ ;phase difference of V2 and V3
  'phase_v2v4c'] ;phase difference of V2 and V4

stop
;--------------------------------
;Non smoothed or detrended versions
tplot,rbspx+['efw_vsvy_v1',$ ;probe potential V1
  'efw_vsvy_v2',$        ;probe potential V2
  'efw_vsvy_v3',$        ;probe potential V3
  'efw_vsvy_v4',$        ;probe potential V4
  'comb12_bp',$  ;smoothed and detrended V1 and V2
  'comb34_bp',$  ;smoothed and detrended V3 and V4
  'comb13_bp',$  ;smoothed and detrended V1 and V3
  'density_smoothed_detrend',$ ;density proxy (V1+V2)/2
  'phase_v1v2c',$ ;phase difference of V1 and V2
  'phase_v3v4c',$ ;phase difference of V3 and V4
  'phase_v1v3c',$ ;phase difference of V1 and V3
  'phase_v2v3c',$ ;phase difference of V2 and V3
  'phase_v2v4c'] ;phase difference of V2 and V4




stop

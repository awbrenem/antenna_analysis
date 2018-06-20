;Plot output .tplot files from efw_probe_potentials_test.pro.
;This is mainly used so that people interested in EFW data quality
;can plot this output.

;Restore .tplot file with the variables
path = '~/Desktop/'
date = '20150315'
probe = 'a'

tplot_restore,filenames='RBSP'+probe+'_'+date+'_probe_potentials_test.tplot'

;Make plots prettier
rbsp_efw_init
rbspx = 'rbsp'+probe+'_'



tplot,['envelope_div_flag',$ ;flag for when E12/E34 ratio varies too far from unity
'e120_flag',$                ;flag indicating probable wake effects
rbspx+'comb12_bp',$  ;smoothed and detrended V1 and V2
rbspx+'comb34_bp',$  ;smoothed and detrended V3 and V4
rbspx+'comb13_bp',$  ;smoothed and detrended V1 and V3
rbspx+'density_smoothed_detrend',$ ;density proxy (V1+V2)/2
rbspx+'phase_v1v2c',$ ;phase difference of V1 and V2
rbspx+'phase_v3v4c',$ ;phase difference of V3 and V4
rbspx+'phase_v1v3c',$ ;phase difference of V1 and V3
rbspx+'phase_v2v3c',$ ;phase difference of V2 and V3
rbspx+'phase_v2v4c'] ;phase difference of V2 and V4

stop
;--------------------------------
;Non smoothed or detrended versions
tplot,[rbspx+'efw_vsvy_v1',$ ;probe potential V1
rbspx+'efw_vsvy_v2',$        ;probe potential V2
rbspx+'efw_vsvy_v3',$        ;probe potential V3
rbspx+'efw_vsvy_v4',$        ;probe potential V4
rbspx+'comb12_bp',$  ;smoothed and detrended V1 and V2
rbspx+'comb34_bp',$  ;smoothed and detrended V3 and V4
rbspx+'comb13_bp',$  ;smoothed and detrended V1 and V3
rbspx+'density_smoothed_detrend',$ ;density proxy (V1+V2)/2
rbspx+'phase_v1v2c',$ ;phase difference of V1 and V2
rbspx+'phase_v3v4c',$ ;phase difference of V3 and V4
rbspx+'phase_v1v3c',$ ;phase difference of V1 and V3
rbspx+'phase_v2v3c',$ ;phase difference of V2 and V3
rbspx+'phase_v2v4c'] ;phase difference of V2 and V4



stop

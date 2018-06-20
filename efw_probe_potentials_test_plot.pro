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

tplot,['envelope_diff_flag',$
'envelope_div_flag',$
'e120_flag','e340_flag',$
'e12_envelopemxmn_test',$
'e34_envelopemxmn_test',$
'e12_e34_envelopemx_comp',$
'e12_e34_envelopemn_comp',$
'envelopemx_diff_comb',$
'envelopemx_div_comb']



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


store_data,['rbsp?_efw_vsvy_v?_smoothed_detrend',$
  'rbsp?_comb??d','euv_comb'],/delete


TPLOT(346):    1 rbspa_efw_vsvy_v1
TPLOT(346):    2 rbspa_efw_vsvy_v2
TPLOT(346):    3 rbspa_efw_vsvy_v3
TPLOT(346):    4 rbspa_efw_vsvy_v4
TPLOT(346):   35 rbspa_density_smoothed_detrend
TPLOT(346):   52 rbspa_phase_v1v2c
TPLOT(396):  36   rbspa_phase_v1v2
TPLOT(346):   57 rbspa_phase_v3v4c
TPLOT(396):  42   line0
TPLOT(396):  43   line90
TPLOT(396):  49   line180red
TPLOT(396):  45   line270
TPLOT(396):  46   line360
TPLOT(396):  41   rbspa_phase_v3v4
TPLOT(346):   53 rbspa_phase_v1v3c
TPLOT(396):  42   line0
TPLOT(396):  48   line90red
TPLOT(396):  44   line180
TPLOT(396):  50   line270red
TPLOT(396):  46   line360
TPLOT(396):  37   rbspa_phase_v1v3
TPLOT(346):   55 rbspa_phase_v2v3c
TPLOT(396):  42   line0
TPLOT(396):  48   line90red
TPLOT(396):  44   line180
TPLOT(396):  50   line270red
TPLOT(396):  46   line360
TPLOT(396):  39   rbspa_phase_v2v3
TPLOT(346):   56 rbspa_phase_v2v4c
TPLOT(396):  42   line0
TPLOT(396):  48   line90red
TPLOT(396):  44   line180
TPLOT(396):  50   line270red
TPLOT(396):  46   line360
TPLOT(396):  40   rbspa_phase_v2v4
TPLOT(346):   63 e120_flag
TPLOT(346):   17 rbspa_comb12_bp
TPLOT(396):  13   rbspa_v1_bp
TPLOT(396):  14   rbspa_v2_bp
TPLOT(346):   22 rbspa_comb34_bp
TPLOT(396):  15   rbspa_v3_bp
TPLOT(396):  16   rbspa_v4_bp
TPLOT(346):   18 rbspa_comb13_bp
TPLOT(396):  13   rbspa_v1_bp
TPLOT(396):  15   rbspa_v3_bp
TPLOT(346):   35 rbspa_density_smoothed_detrend
TPLOT(346):   52 rbspa_phase_v1v2c
TPLOT(396):  36   rbspa_phase_v1v2
TPLOT(346):   57 rbspa_phase_v3v4c
TPLOT(396):  41   rbspa_phase_v3v4
TPLOT(346):   53 rbspa_phase_v1v3c
TPLOT(396):  37   rbspa_phase_v1v3
TPLOT(346):   55 rbspa_phase_v2v3c
TPLOT(396):  39   rbspa_phase_v2v3
TPLOT(346):   56 rbspa_phase_v2v4c
TPLOT(396):  40   rbspa_phase_v2v4

;************************************
23 rbspa_comb12                        rbspa_efw_vsvy_v1 rbspa_efw_vsvy_v2
24 rbspa_comb13                        rbspa_efw_vsvy_v1 rbspa_efw_vsvy_v3
25 rbspa_comb14                        rbspa_efw_vsvy_v1 rbspa_efw_vsvy_v4
26 rbspa_comb23                        rbspa_efw_vsvy_v2 rbspa_efw_vsvy_v3
27 rbspa_comb24                        rbspa_efw_vsvy_v2 rbspa_efw_vsvy_v4
28 rbspa_comb34                        rbspa_efw_vsvy_v3 rbspa_efw_vsvy_v4
29 rbspa_comb12d                       rbspa_efw_vsvy_v1_smoothed_detrend rbspa_efw_vsvy_v2_smoothed_detrend
30 rbspa_comb13d                       rbspa_efw_vsvy_v1_smoothed_detrend rbspa_efw_vsvy_v3_smoothed_detrend
31 rbspa_comb14d                       rbspa_efw_vsvy_v1_smoothed_detrend rbspa_efw_vsvy_v4_smoothed_detrend
32 rbspa_comb23d                       rbspa_efw_vsvy_v2_smoothed_detrend rbspa_efw_vsvy_v3_smoothed_detrend
33 rbspa_comb24d                       rbspa_efw_vsvy_v2_smoothed_detrend rbspa_efw_vsvy_v4_smoothed_detrend
34 rbspa_comb34d                       rbspa_efw_vsvy_v3_smoothed_detrend rbspa_efw_vsvy_v4_smoothed_detrend
35 rbspa_density_smoothed_detrend
36 rbspa_phase_v1v2
37 rbspa_phase_v1v3
38 rbspa_phase_v1v4
39 rbspa_phase_v2v3
40 rbspa_phase_v2v4
41 rbspa_phase_v3v4
42 line0
43 line90
44 line180
45 line270
46 line360
47 line0red
48 line90red
49 line180red
50 line270red
51 line360red
52 rbspa_phase_v1v2c                   line0 line90 line180red line270 line360 rbspa_phase_v1v2
53 rbspa_phase_v1v3c                   line0 line90red line180 line270red line360 rbspa_phase_v1v3
54 rbspa_phase_v1v4c                   line0 line90red line180 line270red line360 rbspa_phase_v1v4
55 rbspa_phase_v2v3c                   line0 line90red line180 line270red line360 rbspa_phase_v2v3
56 rbspa_phase_v2v4c                   line0 line90red line180 line270red line360 rbspa_phase_v2v4
57 rbspa_phase_v3v4c                   line0 line90 line180red line270 line360 rbspa_phase_v3v4
58 e12
59 e34
60 e120
61 e120_mat
62 e120_amp
63 e120_flag
64 e340
65 e340_mat
66 e340_amp
67 e340_flag
68 e12_smoothed
69 e12_detrend
70 e34_smoothed
71 e34_detrend
72 euv_comb                            e12_detrend e34_detrend
73 e12_envelopemx
74 e12_envelopemn
75 e34_envelopemx
76 e34_envelopemn
77 envelope_diff_flag
78 envelope_div_flag
79 e12_envelopemxmn_test               e12_envelopemx e12_envelopemn e12_detrend
80 e34_envelopemxmn_test               e34_envelopemx e34_envelopemn e34_detrend
81 e12_e34_envelopemx_comp             e12_envelopemx e34_envelopemx
82 e12_e34_envelopemn_comp             e12_envelopemn e34_envelopemn
83 e12_e34_envelopemx_diff
84 e12_e34_envelopemn_diff
85 e12_e34_envelopemx_div
86 e12_e34_envelopemn_div
87 diffline
88 divline
89 envelopemx_diff_comb                e12_e34_envelopemx_diff e12_e34_envelopemn_diff diffline
90 envelopemx_div_comb                 e12_e34_envelopemx_div e12_e34_envelopemn_div divline


stop

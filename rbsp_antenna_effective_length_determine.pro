;Crib sheet for finding the RBSP antenna effective length. 
;There's a pretty systematic error in the vxB subtraction results. 
;I'm finding what this error is by finding the slope of the plot of
;actual Esvy mgse data vs the vxb removed field. 
;I then "fix" the original Esvy MGSE data with the equations below. Finally, I
;do a new vxB subtraction. This results in an improvement over the original 
;vxB subtracted field. 

;The Esvy data are adjusted by the effective antenna length through this equation:
;Esvy_adj_y = Esvy_y + Esvy_y*(1-slope_y)
;Esvy_adj_z = Esvy_z + Esvy_z*(1-slope_z)


;Created by Aaron W Breneman
;University of Minnesota  2013-01-24
;awbrenem@gmail.com



;Good days to test

;- [ ] 2012-10-20 - very quiet - eclipse day
;- [ ] 2012-12-14 - very quiet - no eclipse
;- [ ] 2012-12-22 - very quiet - no eclipse
;- [ ] 2012-12-27
;- [ ] 2013-01-02 - possible Msheath
;- [ ] 2013-02-01 - possible Msheath




;Results from different days:

;------------------------------------
;2012-10-13
;RBSP-A
;t0 = time_double('2012-10-13/16:00')
;t1 = time_double('2012-10-13/24:00')
;date = '2012-10-13'
;timespan,date
;probe = 'a'
;slopey = 0.94154335
;intercepty = 0.44745866
;slopez = 0.94483345
;interceptz = 0.67359717
;------------------------------------
;2012-10-19 (full-day)
;RBSP-A,B
;t0 = time_double('2012-10-19/00:00')
;t1 = time_double('2012-10-19/24:00')

;RBSP-A
;slopey = 0.94137313
;intercepty = 0.49965222
;slopez = 0.93839995
;interceptz = 0.55242079

;RBSP-B
;slopey = 0.94515762
;intercepty = 0.56782032
;slopez = 0.94470268
;interceptz = 0.47807721

;-----------------------------------

;Very quiet day - eclipse

;date = '2012-10-20'
;probe = 'a'
;rbspx = 'rbsp'+probe
;duration = 1	; days.
;timespan, date, duration

;;Times for low density Msphere
;t0 = time_double(date+'/09:30')
;t1 = time_double(date+'/13:00')

;;Times for V b/t 0 and -7 (plasmasphere)
;t0 = time_double(date+'/06:20')
;t1 = time_double(date+'/06:35')


;FULL DAY FIT VALUES (minus eclipse)
;RBSP-A
;slopey =       0.94553640
;intercepty =       0.36063051
;slopez =       0.93540665
;interceptz =       0.73798039

;RBSP-B
;slopey =       0.94812763
;intercepty =       0.43458006
;Not completely linear for z
;-------------------------------------
;date = '2012-10-21'
;probe = 'a'
;t0 = time_double(date+'/02:00')
;t1 = time_double(date+'/09:32')

;probe = 'b'
;t0 = time_double(date+'/02:00')
;t1 = time_double(date+'/07:05')

;RBSP-A  
;slopey =       0.94620636
;intercepty =       0.42738234
;slopez =       0.94985750
;interceptz =       0.61703434
  
;RBSP-B
;slopey =       0.95016421
;intercepty =       0.34081106
;slopez =       0.94405647
;interceptz =       0.92197044
;---------------------------------------
;date = '2012-10-22'
;probe = 'a'
;t0 = time_double(date+'/05:00')
;t1 = time_double(date+'/12:28')
	
;slopey =       0.94670716
;intercepty =       0.39294199
;slopez =       0.96193262
;interceptz =       0.47578776

            
;------------------------------------
;2012-11-01
;RBSP-A
;t0 = time_double('2012-11-01/00:00')
;t1 = time_double('2012-11-01/24:00')

;slopey = 0.93537717
;intercepty = 2.0413767
;slopez = 0.93913572
;interceptz = 0.91490549

;t0 = time_double('2012-11-01/00:00')
;t1 = time_double('2012-11-01/07:00')

;slopey = 0.94855319
;intercepty = 0.53938971
;slopez = 0.93464284
;interceptz =  0.65869053
;-------------------------------------

date = '2012-12-23'
probe = 'a'
;t0 = time_double(date + '/15:20')
;t1 = time_double(date + '/16:00')
t0 = time_double(date + '/05:00')
t1 = time_double(date + '/07:00')

;Values without removal of corotation field.
;slopey = 0.95202853
;intercepty = 1.07
;slopez = 0.97336026
;interceptz = 0.307

;Values with removal of corotation field
;slopey = 1.0138
;intercepty = -0.05877
;slopez = 1.03269
;interceptz = 0.0904




;Try the individual perigee passes separately

;Very quiet days to test
;1025, 1026, 1029, 1030, 1031, 1103

timespan,date
rbspx = 'rbsp' + probe


rbsp_efw_init

v = -1*indgen(200)/(10.)


;potential to density
;For October and most of November
N=7927.28*exp(3.816*v)+232.55*exp(0.922*v)+5.57
;The density fit for Nov 27 and after is
N=5795.6*exp(4.304*v)+165.01*exp(0.926*v)+2.29
;where v is the negative of sc potential, and the density, N, is in cm^-3

plot,v,N,xrange=[-20,0],/ylog,ytitle='Density',xtitle='inverse SC pot'


;Looks like Scott's density fit doesn't work for the lowest sc potentials. 


;Region		density		temp (eV)
;Lobes 		0.001-0.1	100-500
;Psheet		0.1-1		200-1000
;MS (outer)	0.1-10		1-100
;Msheath	10-100		100
;Psphere	100-1000	1			(v dominated by sc motion)



;date = '2012-10-22'
;probe = 'a'
;rbspx = 'rbsp'+probe
;duration = 1	; days.
;timespan, date, duration
;;Grab a subset of the data to fit line to.
;t0 = time_double(date+'/05:00')
;t1 = time_double(date+'/12:28')

rbsp_load_spice_kernels


;Determine corotation Efield
rbsp_corotation_efield,probe,date,no_spice_load=no_spice_load


rbsp_load_eclipse_predict,probe,date	
;get_data,'eclipse_umbra_start',data=eus
;get_data,'eclipse_umbra_stop',data=eue
;get_data,'eclipse_penumbra_start',data=eps
;get_data,'eclipse_penumbra_stop',data=epe
get_data,rbspx+'_umbra',data=eu
get_data,rbspx+'_peumbra',data=ep
	
;eclipsetimes = date + '/' + ['03:31:40.781','03:31:50.181','03:31:50.181','04:07:55.083',$
;	'04:07:55.083','04:08:16.389','12:30:29.213','12:30:38.661','12:30:38.661','13:06:46.005',$
;'13:06:46.005','13:07:07.397','21:29:17.721','21:29:27.217','21:29:27.217','22:05:36.930',$
;'22:05:36.930','22:05:58.409']

	

tplot_options,'xmargin',[20.,15.]
tplot_options,'ymargin',[3,6]
tplot_options,'xticklen',0.08
tplot_options,'yticklen',0.02
tplot_options,'xthick',2
tplot_options,'ythick',2
tplot_options,'labflag',-1	

 


rbsp_load_emfisis,probe=probe,coord='gse',cadence='1sec',level='l3'  ;load this for the mag model subtract
get_data,rbspx+'_emfisis_l3_1sec_gse_Mag',data=dat,dlimits=dlim,limits=lim	
store_data,rbspx+'_mag_gse',data=dat,dlimits=dlim,limits=lim




rbsp_load_efw_waveform,probe=probe,type='calibrated',datatype='vsvy'
spinperiod = 10.8
rbsp_downsample,rbspx+'_efw_vsvy',1/spinperiod,/nochange	

;Create the density proxy variable
split_vec,rbspx+'_efw_vsvy'

get_data,rbspx+'_efw_vsvy_0',data=v1
get_data,rbspx+'_efw_vsvy_1',data=v2
get_data,rbspx+'_efw_vsvy_2',data=v3
get_data,rbspx+'_efw_vsvy_3',data=v4
get_data,rbspx+'_efw_vsvy_4',data=v5
get_data,rbspx+'_efw_vsvy_5',data=v6

dproxy = (v1.y + v2.y + v3.y + v4.y)/4.
dproxy12 = (v1.y + v2.y)/2.
dproxy34 = (v3.y + v4.y)/2.
dproxy56 = (v5.y + v6.y)/2.

store_data,'density_proxy',data={x:v1.x,y:dproxy}
store_data,'density_proxy12',data={x:v1.x,y:dproxy12}
store_data,'density_proxy34',data={x:v1.x,y:dproxy34}
store_data,'density_proxy56',data={x:v1.x,y:dproxy56}

store_data,'dp_comb',data=['density_proxy','density_proxy12','density_proxy34','density_proxy56']
options,'dp_comb','colors',[0,1,2,3]
tplot,['dp_comb',rbspx+'_efw_vsvy_0']
if is_struct(eu) then timebar,eu.x,color=50
if is_struct(eu) then timebar,eu.x+eu.y,color=50
if is_struct(ep) then timebar,ep.x,color=50
if is_struct(ep) then timebar,ep.x+ep.y,color=50




rbsp_load_efw_esvy_mgse,probe=probe,/no_spice_load

;Remove corotation Efield
dif_data,rbspx+'_efw_esvy_mgse',rbspx+'_E_coro_mgse',newname=rbspx+'_efw_esvy_mgse'


split_vec,rbspx+'_efw_esvy_mgse'

tplot,[rbspx+'_efw_esvy_mgse_x',rbspx+'_efw_esvy_mgse_y',rbspx+'_efw_esvy_mgse_z']




;Get antenna pointing direction and stuff
rbsp_load_state,probe=probe,/no_spice_load,datatype=['spinper','spinphase','mat_dsc','Lvec'] 

rbsp_efw_position_velocity_crib,/no_spice_load,/noplot
get_data,rbspx+'_spinaxis_direction_gse',data=wsc_GSE	


;rbsp_load_spice_state,probe=probe,coord='gse',/no_spice_load  
get_data,rbspx+'_state_pos_gse',data=pos_gse


;Create velocity magnitude variables
get_data,rbspx+'_state_vel_gse',data=vel
	

cotrans,rbspx+'_state_pos_gse',rbspx+'_state_pos_gsm',/GSE2GSM	





;time2=time_double(date) ; first get unix time double for beginning of day
;
;;Grab first and last time of day
;time2 = [time2,time2+86399.]
;time3=time_string(time2, prec=6) ; turn it back into a string for ISO conversion
;strput,time3,'T',10 ; convert TPLOT time string 'yyyy-mm-dd/hh:mm:ss.msec' to ISO 'yyyy-mm-ddThh:mm:ss.msec'
;cspice_str2et,time3,et2 ; convert ISO time string to SPICE ET
;
;
;cspice_pxform,'RBSP'+strupcase(probe)+'_SCIENCE','GSE',et2[0],pxform1
;cspice_pxform,'RBSP'+strupcase(probe)+'_SCIENCE','GSE',et2[1],pxform2
;
;wsc1=dblarr(3)
;wsc1[2]=1d
;wsc_GSE1=dblarr(3)
;
;wsc2=dblarr(3)
;wsc2[2]=1d
;wsc_GSE2=dblarr(3)
;
;; Calculate the modified MGSE directions	
;wsc_GSE1 = pxform1 ## wsc1  ;start of day
;wsc_GSE2 = pxform2 ## wsc2  ;end of day	
;
;;Puts velocity and mag data in MGSE
;rbsp_gse2mgse,rbspx+'_state_vel_gse',reform(wsc_GSE1),newname=rbspx+'_vel_mgse'

get_data,rbspx+'_mag_gse',data=tmpp2
wsc_GSE_tmp2 = [[interpol(wsc_GSE[0,*],time_double(time3),tmpp2.x)],$
			   [interpol(wsc_GSE[1,*],time_double(time3),tmpp2.x)],$
			   [interpol(wsc_GSE[2,*],time_double(time3),tmpp2.x)]]

rbsp_gse2mgse,rbspx+'_mag_gse',reform(wsc_GSE_tmp2),newname=rbspx+'_mag_mgse'

options,rbspx+'_mag_mgse','colors',[2,4,6]
	


;Remove vxB field
rbsp_vxb_subtract,rbspx+'_vel_mgse',rbspx+'_mag_mgse',rbspx+'_efw_esvy_mgse'
get_data,'Esvy_mgse_vxb_removed',data=dd
store_data,'E_vxb_removed_orig',data=dd
store_data,'Esvy_mgse_vxb_removed',/delete
split_vec,'E_vxb_removed_orig'

ylim,['E_vxb_removed_orig_x','E_vxb_removed_orig_y','E_vxb_removed_orig_z'],-20,20
tplot,['density_proxy','E_vxb_removed_orig_x','E_vxb_removed_orig_y','E_vxb_removed_orig_z']
if is_struct(eu) then timebar,eu.x,color=50
if is_struct(eu) then timebar,eu.x+eu.y,color=50
if is_struct(ep) then timebar,ep.x,color=50
if is_struct(ep) then timebar,ep.x+ep.y,color=50






ylim,['E_vxb_removed_orig_x','E_vxb_removed_orig_y','E_vxb_removed_orig_z'],-20,20
tplot,['density_proxy','E_vxb_removed_orig_x','E_vxb_removed_orig_y','E_vxb_removed_orig_z']
if is_struct(eu) then timebar,eu.x,color=50
if is_struct(eu) then timebar,eu.x+eu.y,color=50
if is_struct(ep) then timebar,ep.x,color=50
if is_struct(ep) then timebar,ep.x+ep.y,color=50

	


ylim,[rbspx+'_efw_esvy_mgse_y','vxb_y'],-300,300
ylim,'E_vxb_removed_orig_y',-20,20
ylim,[rbspx+'_efw_esvy_mgse_z','vxb_z'],-300,300
ylim,'E_vxb_removed_orig_z',-20,20
tplot,['E_vxb_removed_orig_y',rbspx+'_efw_esvy_mgse_y','vxb_y']
tplot,['E_vxb_removed_orig_z',rbspx+'_efw_esvy_mgse_z','vxb_z']
if is_struct(eu) then timebar,eu.x,color=50
if is_struct(eu) then timebar,eu.x+eu.y,color=50
if is_struct(ep) then timebar,ep.x,color=50
if is_struct(ep) then timebar,ep.x+ep.y,color=50




;*****************
;TEST: ADD IN EFFECTIVE LENGTH AND COMPARE

	;Apply crude antenna effective length correction to y and z MGSE values
	get_data,'E_vxb_removed_orig',data=dat,dlimits=dlim,limits=lim
	dat.y[*,1] = dat.y[*,1] + dat.y[*,1]*(1-1.014)
	dat.y[*,2] = dat.y[*,2] + dat.y[*,2]*(1-1.033)
;	dat.y[*,1] = dat.y[*,1] + dat.y[*,1]*(1-1.9)
;	dat.y[*,2] = dat.y[*,2] + dat.y[*,2]*(1-1.9)
	store_data,'E_vxb_removed_orig22',data=dat,dlimits=dlim,limits=lim

	split_vec,'E_vxb_removed_orig22'
ylim,['E_vxb_removed_orig_y','E_vxb_removed_orig22_y',$
	'E_vxb_removed_orig_z','E_vxb_removed_orig22_z'],-5,5

tplot,['E_vxb_removed_orig_y','E_vxb_removed_orig22_y',$
	'E_vxb_removed_orig_z','E_vxb_removed_orig22_z']


;Values with removal of corotation field
;slopey = 1.0138
;intercepty = -0.05877
;slopez = 1.03269
;interceptz = 0.0904




;******************










;Choose data points by trange

	get_data,rbspx+'_efw_esvy_mgse_y',data=dy
	get_data,'vxb_y',data=my

	get_data,rbspx+'_efw_esvy_mgse_z',data=dz
	get_data,'vxb_z',data=mz

	gooy = where((dy.x ge t0) and (dy.x le t1))
	gooz = where((dz.x ge t0) and (dz.x le t1))




;Fit a straight line to the real (Esurvey) data vs model (vxB) field data
;mey = replicate(8.0,n_elements(gooy))
;mez = replicate(8.0,n_elements(gooz))

resulty = linfit(my.y[gooy],dy.y[gooy]);,chisqr=chsqy,measure_errors=mey)
resultz = linfit(mz.y[gooz],dz.y[gooz]);,chisqr=chsqz,measure_errors=mez)


;define independent variables for plotting
xvalsy = dindgen(n_elements(gooy))*max(my.y[gooy])/(n_elements(gooy)-1)
xvalsz = dindgen(n_elements(gooz))*max(mz.y[gooz])/(n_elements(gooz)-1)

;Equation of best-fit line
fitliney = xvalsy*resulty[1] + resulty[0]
fitlinez = xvalsz*resultz[1] + resultz[0]

;Line with slope=1
idealliney = xvalsy*1 + 0
ideallinez = xvalsz*1 + 0


popen,'~/Desktop/eff_antenna_length_fit_RBSP'+probe+'_'+date

	xr = [-20,160]  ;good for perigee
;	xr = [-0.5,0.5]

	!p.charsize = 2
	!p.multi = [0,0,3]

	plot,my.y[gooy],dy.y[gooy],psym=4,yrange=xr,xrange=xr,ystyle=1,xstyle=1,$
		xtitle='vxB MGSEy (mV/m)',ytitle='Esvy MGSEy (mV/m)',title=strupcase(rbspx) + ' ' + date
	oplot,xvalsy,fitliney,color=250,thick=2
	oplot,xvalsy,idealliney,color=50,thick=2

	plot,mz.y[gooy],dz.y[gooy],psym=3,yrange=xr,xrange=xr,ystyle=1,xstyle=1,$
		xtitle='vxB MGSEz (mV/m)',ytitle='Esvy MGSEz (mV/m)',title=strupcase(rbspx) + ' ' + date
	oplot,xvalsz,fitlinez,color=250,thick=2
	oplot,xvalsz,ideallinez,color=50,thick=2

	plot,my.y[gooy],dy.y[gooy],psym=3,yrange=xr,xrange=xr,ystyle=1,xstyle=1,$
		xtitle='vxB MGSEy,z (mV/m)',ytitle='Esvy MGSEy,z (mV/m)',title=strupcase(rbspx) + ' ' + date,/nodata
	oplot,xvalsy,fitliney,color=210,thick=2
	oplot,xvalsz,fitlinez,color=250,thick=2
	oplot,xvalsz,ideallinez,color=50
 

	xyouts,0.12,0.95,'Blue line is slope of 1',/norm,charsize=0.8
	xyouts,0.12,0.92,'Red and orange lines are best fit lines for y and z MGSE',/norm,charsize=0.8
	xyouts,0.12,0.62,'Ey MGSE slope = '+strtrim(resulty[1],2),/norm,charsize=0.8
	xyouts,0.12,0.57,'Ey MGSE intercept = '+strtrim(resulty[0],2),/norm,charsize=0.8
	xyouts,0.12,0.52,'Ez MGSE slope = '+strtrim(resultz[1],2),/norm,charsize=0.8
	xyouts,0.12,0.47,'Ez MGSE intercept = '+strtrim(resultz[0],2),/norm,charsize=0.8
	xyouts,0.12,0.28,time_string(t0)+' to '+time_string(t1),/norm,charsize=0.8


pclose


print,'slopey = ',resulty[1]
print,'intercepty = ',resulty[0]
print,'slopez = ',resultz[1]
print,'interceptz = ',resultz[0]




get_data,rbspx+'_efw_esvy_mgse_x',data=dx



;Multiply the real EFW data by the slope. This may be the antenna effective length correction.
dy_adj = dy.y + dy.y*(1-resulty[1])
dz_adj = dz.y + dz.y*(1-resultz[1])

;Here's the amount we add to the Esvy data to account for the new effective antenna length
store_data,'dy_adj',data={x:dy.x,y:dy.y*(1-resulty[1])}
store_data,'dz_adj',data={x:dz.x,y:dz.y*(1-resultz[1])}
tplot_options,'title','Amount we add to the Esvy data!Cto account for the new determined antenna length'
options,['dy_adj','dz_adj'],'ysubtitle','[mV/m]'
options,'dy_adj','ytitle','delta-y MGSE'
options,'dz_adj','ytitle','delta-z MGSE'
tplot,['dy_adj','dz_adj']


;;Esvy data corrected by the new effective antenna length.
dnew = [[dx.y],[dy_adj],[dz_adj]]
store_data,'Emgse_mgse_adj',data={x:dy.x,y:dnew}
;tplot,[rbspx+'_efw_esvy_mgse','Emgse_mgse_adj']


;Now calculate the vxB removal using the "corrected" Esvy data

	rbsp_vxb_subtract,rbspx+'_vel_mgse',rbspx+'_mag_mgse','Emgse_mgse_adj'
	get_data,'Esvy_mgse_vxb_removed',data=dd
	store_data,'E_vxb_removed_adj',data=dd
	split_vec,'E_vxb_removed_adj'
	ylim,['E_vxb_removed_orig_y','E_vxb_removed_adj_y'],-20,20
	tplot_options,'title','Esvy MGSE y and z data before (black) and after!Ceffective antenna length correction'
	options,['E_vxb_removed_orig_y','E_vxb_removed_adj_y'],'ysubtitle','[mV/m]'
	options,'E_vxb_removed_orig_y','ytitle','Esvy MGSEy!CvxB removed!Coriginal'
	options,'E_vxb_removed_adj_y','ytitle','Esvy MGSEy!CvxB removed!Ccorrected'

	ylim,['E_vxb_removed_orig_z','E_vxb_removed_adj_z'],-20,20
;	tplot_options,'title','Esvy MGSE z data before and after!Ceffective antenna length correction'
	options,['E_vxb_removed_orig_z','E_vxb_removed_adj_z'],'ysubtitle','[mV/m]'
	options,'E_vxb_removed_orig_z','ytitle','Esvy MGSEz!CvxB removed!Coriginal'
	options,'E_vxb_removed_adj_z','ytitle','Esvy MGSEz!CvxB removed!Ccorrected'



rbsp_downsample,['E_vxb_removed_orig_y','E_vxb_removed_adj_y',$
				 'E_vxb_removed_orig_z','E_vxb_removed_adj_z'],1/11.4
store_data,'Ey_comp_DS',data=['E_vxb_removed_orig_y_DS','E_vxb_removed_adj_y_DS']
store_data,'Ez_comp_DS',data=['E_vxb_removed_orig_z_DS','E_vxb_removed_adj_z_DS']
options,'Ey_comp_DS',colors=[0,250]
options,'Ez_comp_DS',colors=[0,50]
options,'Ey_comp_DS','ytitle','Esvy MGSEy!CvxB removed'
options,'Ez_comp_DS','ytitle','Esvy MGSEz!CvxB removed'
options,'Ey_comp_DS','ysubtitle','[mV/m]'
options,'Ez_comp_DS','ysubtitle','[mV/m]'



popen,'~/Desktop/eff_antenna_length_fit_RBSP'+probe+'_'+date+'_comparison',/landscape
	!p.charsize = 0.9
	tplot,['Ey_comp_DS','Ez_comp_DS']
	timebar,eclipsetimes,color=210
pclose


;tplot,['E_vxb_removed_orig_z','E_vxb_removed_adj_z']




store_data,'combined',data=['vxb_y',rbspx+'_efw_esvy_mgse_y']

;tplot,['E_vxb_removed_orig_y','E_vxb_removed_adj_y','vxb_y',rbspx+'_efw_esvy_mgse_y']
tplot,['E_vxb_removed_orig_y','E_vxb_removed_adj_y','combined']




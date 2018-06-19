;Calculate the antenna response for electric field antennas in space plasmas

;PATDSGAIN computes:
;	THIS ROUTINE COMPUTES GAINS OF THE TDS AND PREAMP
;	THE INPUT IS TAKEN FROM THE PREAMP INPUT,
;	ANTENNA IMPEDANCE NOT TAKEN INTO ACCOUNT

;Written by Aaron W Breneman, 2017


;dens = density in cm-3
;temp = temp in eV



;------------------------------------
;TYPICAL SOLAR WIND PARAMETERS AT 1AU
;------------------------------------

;fast SW
;Tp = 10-20 eV
;Te = 5-20 eV
;n = 1-10 cm-3

;slow SW
;Tp = 5-20 eV
;Te = 5-20 eV
;n = 5-25 cm-3

pro spacecraft_sheath_values,dens,temp,eclipse=eclipse,cylindrical=cyl


;--------------
;Constants
;--------------
epo    = 8.854187817d-12   ; -Permittivity of free space (F/m)
e     = 1.60217733d-19    ; -Fundamental charge (C)
kB     = 1.380658d-23      ; -Boltzmann Constant (J/K)
c      = 2.99792458d8      ; -Speed of light in vacuum (m/s)
me_ev  = 0.5109906d6       ; -Electron mass in eV/c^2
mp_ev  = 938.27231d6       ; -Proton mass in eV/c^2
mi = 938d6/(c^2)  ;eV/c^2
me  = 0.51d6/(c^2)       ; -Electron mass in eV/c^2




if debye lt a then begin
	print,'*********WARNING***************************************************************'
	print,"The Debye length is smaller than the antenna radius...Approximations don't hold"
	print,"Manually setting Debye = antenna radius"
	print,'*******************************************************************************'
	debye = a
endif


if KEYWORD_SET(stereo) then begin
	ibias = 0.
	As = 0.45d ;m^2 Antenna surface area from Bale08
	L = 6d  ;m (antenna length)
	a = 0.023d      ;m  (average antenna radius from Bale08)
	r = 0.      ;probe radius
endif




;-------------------------------
;SHEATH RESISTANCE (DAYSIDE)
;-------------------------------
;From Gurnett, Eqn 6. principles of plasma wave instrument design
;Rs = Up/(Ie - Ii - Ibias)
;Ibias = 0
;Up ~ 1.5 volt. (See Cauffman and Gurnett72, but I think this value is pretty static)
;Ii = 0  (electron current should dominate)
;Ie -> e- thermal current (See Kersten eqn 48, Kellogg09 eqn 2 or Mott and Langmuir eqn 24,25)
;This eqn is confirmed by Stuart's email
;If it was sunlit, then a good estimate of the sheath resistance is $R_s=kT_p/I_e$ where $T_p$
;is the photoelectron energy (approx 1V) and $I_e$ is the thermal electron current (see Gurnett eqn 6).



;------------------
;SHEATH CAPACITANCE
;------------------
;From Gurnett, Eqn 7. principles of plasma wave instrument design
;This assumes that Debye >> a (antenna radius). This may not be true

;Cs = 2*!pi*epo*(L/2)/(ln(Debye*a)) --> units of Farads = Ampere*s/Volts



;-----------------------------------------------------------------
;STEREO G(w) for the antenna interaction with the plasma sheath.
;See unshaded part of Bale08 figure 2
;This equation is completely general
;Note: The 32pF value for the base capacitance is the value that makes
;	STEREO antennas act like dipoles (monopole to infinite grounded plane)
;	This value is also taken into account with the PATDS program, so
;-----------------------------------------------------------------

;To transform from ANTENNA voltages to electric fields we use
;the effective lengths and angles from Bale08 (just under Fig 14).
;The question for low freq waves is, what are the effective lengths
;and angles? These were determined for high freqs in Bale08 but
;not for low freqs. They were determined from the wire Gratz test
;using base caps with 90 pF.



;s = sheath
;b = base (voltage divider effect that occurs b/c of the way the antennas are attached to sc)

;Rs = 0.75d   ;Mohms
;Cs = 40d 	;pF

;Rb = 1d		;Mohms
;Cb = 67d		;pF  (actually stray capacitance from Bale08
;				;Cstray = Cb + Ccoax_cables + Cpreamp_input
;				;Cstray = 32 + 13 + 22 = 67
;				;The values of the coax and preamp input should be stable from low
;				;freqs up to GHz freqs.



;------------------
;DEFINE SOME VALUES
;------------------



;-------------------------------
;SHEATH RESISTANCE (NIGHTSIDE)
;-------------------------------
;From Gurnett, Eqn 5. principles of plasma wave instrument design
;for the nightside PS when the sc is in shadow and is charged negative.
;Rs = Ue/(Ip + Ii + Ibias)  -> units of (V*s/C=V/Amp=Ohms)

;In this case Ip=0 (no photoionization current)
;Ibias = 0
;Ue = (kTe/e)
;Ii = Isat = Ae*n_o/2 * sqrt(kTe/mi)  ;Let the ion saturation current be the Bohm current
;
;This equation is confirmed by an email from Stuart Bale:
;If it's dark, then it's more like $R_s = kT_e/I_i$  - the thermal electron energy over the ion thermal current.


;ION CURRENTS TEND TO BE ONLY FEW % OF E- CURRENTS EXCEPT DURING TIMES OF EXTREME
;PLASMA FLOWS IN SOLAR WIND OR IN THE PLASMASHEET (PEDERSON08 [8])

;IBIAS TENDS TO DOMINATE Ie IN LOW DENSITY (LOBE) CONDITIONS. OTHERWISE THE TWO
;CAN BE COMPARABLE.


;--------------------------------------------------------
;GENERAL QUANTITIES (independent of sc or antenna type)
;--------------------------------------------------------

no = dens * (100d^3)
kTe = temp
debye = 7.43d2*sqrt(kTe)/sqrt(no/(100^3))/100.			  ;Debye length (m)


;Electron thermal current
;--Pedersen08
C = 2.68d-14  ;A*m^3/sqrt(V)
Vp_np = (Vp - Vnp) ;probe potential relative to local plasma
Ve = kTe/e  ;e- thermal velocity
Iep0 = C*no*sqrt(Ve)*Ap
Ie = Ie0*(1 + Vp_np/Ve)










;------------------------------------------------------------------------------

;Cs is same for eclipse and sunlit
if type eq 'cylindrical' then begin

	if keyword_set(eclipse) then begin   ;SHADOWED
;		Rs = Ue/(Ie + Ii + Ibias)   ;general eqn

		Ie = 0.  ;(See Kersten eqn 48, Kellogg09 eqn 2 or Mott and Langmuir eqn 24,25)
		Ii = (As*e*no/2.)*sqrt(kTe/mi)  ;IS THIS ONLY FOR STEREO?
		if keyword_set(stereo) then Rs = (kTe)/Ii
		Cs = 2*!pi*epo*(L/2.)/alog(debye/a) ;general eqn

	endif else begin  ;SUNLIT
	 ;Rs = Up/(Ie - Ii - Ibias)  ;general eqn

		Ii = 0.
		Ie = (As*e*no/2.)*sqrt(kTe/me)   ;IS THIS ONLY FOR STEREO???
		Up = 1.5 ;static
		Rs = Up/Ie
		Cs = 2*!pi*epo*(L/2.)/alog(Debye/a) ;general eqn
	endif
endif
;Ie -> e- thermal current (See Kersten eqn 48, Kellogg09 eqn 2 or Mott and Langmuir eqn 24,25)
;This eqn is confirmed by Stuart's email
;If it was sunlit, then a good estimate of the sheath resistance is $R_s=kT_p/I_e$ where $T_p$
;is the photoelectron energy (approx 1V) and $I_e$ is the thermal electron current (see Gurnett eqn 6).


if type eq 'doubleprobe' then begin

;General eqn for Ie to spherical probe: Pederson08 eqn 2



	if keyword_set(eclipse) then begin  ;SHADOWED
		;Rs = Ue/(Iph + Ii + Ibias)   ;general eqn for ion sheath
		Ue = kTe/e
		Ii = 0.  ;I DON'T THINK THIS IS TRUE FOR SHADOWED CONDITIONS
		Rs = Ue/(Ie + Ii + Ibias)
		Cs = 4*!pi*epo*r*(1 + r/debye) ;general eqn

	endif else begin  ;SUNLIT
		;Rs = Up/(Ie - Ii - Ibias)  ;general eqn for e- sheath

		Up = 1.5  ;volts - characteristic energy of photoelectron spectrum (Cauffmann and Gurnett72)
		Ii = 0.
		Ie = -4*!pi*(r^2)*no*e*sqrt(kTe/(2*!pi*me))   ;Mozer16

		;Or, we can write Ie in terms of potential as (Langmuir probe eqn)
		Ies = e*no*vthe*A/4.
		Ie = Ies*exp(-1*(Vp-Vb)/kTe)


		Rs = Up/(Ie - Ii - Ibias)
		Cs = 4*!pi*epo*r*(1 + r/debye) ;general eqn

	endif


endif



print,'-------------------------------------------------------------------'
print,'Rs = ',Rs/1d6,' MegaOhms'
print,'Cs = ',cs*1d12,' pF'
print,'Debye length (m) = ',debye
print,'-------------------------------------------------------------------'


end

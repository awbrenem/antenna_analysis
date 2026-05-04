#Langmuir probe analysis

"""
Determine currents collected by a Langmuir probe in various plasma regimes (unmagnetized, magnetized, thermal, nonthermal, mesothermal) using the equations from Barjatya+2007 thesis.

Useful reference is Barjatya+2007 thesis
"""

import numpy as np


#------------------------------------------------------
#Define constants
#------------------------------------------------------

k = 1.38e-23 #Boltzmann constant (m^2 kg s^-2 K^-1)
q = 1.6e-19 #elementary charge (C)
me = 9.11e-31 #electron mass (kg)
mi = 1.67e-27 #ion mass (kg)
eo = 8.85e-12 #permittivity of free space (s^4 A^2 m^-3 kg^-1)
Vp = 0 #plasma potential (V)
A = 1e-4 #probe surface area (cm-2)


#VERIFIED
def debye_length(ne, Te):
    TeK = Te * 11604.5 #convert from eV to K
    ne = ne * 1e6 #convert from cm^-3 to m^-3
    #return np.sqrt(eo * k * TeK / (2 * np.pi * me)) 
    return np.sqrt(eo * k * TeK / (ne * q**2)) #m


def currentPhotoelectron(Iph0, Te, V, Vp, r):
    #Photoelectron current
    #r is the probe radius (cm)
    #Iph0 is the photoelectron current density at the probe surface (nA/cm^2). 
    #   This depends on the material properties of the probe and the solar UV flux, 
    #   and can be estimated from models or measurements (e.g. Fahleson+74).

    TeK = Te * 11604.5 #convert from eV to K
    return np.pi * r**2 * Iph0 * np.exp(-q * (V - Vp) / (k * TeK)) #nA


#------------------------------------------------------
#Non-drifting, unmagnetized, Maxwellian plasmas, collisionless
#------------------------------------------------------

def currentThermalElectron(ne, Te, A):
    #[Barjatya +2007 thesis Eqn 1.2]
    A = A * 1e-4 #convert from cm^2 to m^2
    ne = ne * 1e6 #convert from cm^-3 to m^-3
    TeK = Te * 11604.5 #convert from eV to K
    return 1e9 * ne * q * A * np.sqrt((k * TeK) / (2 * np.pi * me)) #nA 

def currentThermalIon(ni, Ti, A):
    #[Barjatya +2007 thesis Eqn 1.2]
    A = A * 1e-4 #convert from cm^2 to m^2
    ni = ni * 1e6 #convert from cm^-3 to m^-3
    TiK = Ti * 11604.5 #convert from eV to K
    return 1e9 * ni * q * A * np.sqrt((k * TiK) / (2 * np.pi * mi)) #nA



def currentElectronRetardation(ne, Te, V, Vp, A):
    #[Barjatya +2007 thesis Eqn 1.4]
    #Electron retardation region current (NOTE: e- current dominates ion current here)
    A = A * 1e-4 #convert from cm^2 to m^2
    ne = ne * 1e6 #convert from cm^-3 to m^-3
    TeK = Te * 11604.5 #convert from eV to K
    Ithe = currentThermalElectron(ne, Te, A)
    return Ithe * np.exp(q * (V-Vp) / (k * TeK))  #nA







def currentElectronSaturation(ne, Te, V, Vp, A, b):
    #[Barjatya +2007 thesis Eqn 1.5]
    #Electron saturation current (NOTE: e- current dominates ion current here)
    A = A * 1e-4 #convert from cm^2 to m^2
    ne = ne * 1e6 #convert from cm^-3 to m^-3
    TeK = Te * 11604.5 #convert from eV to K
    Ithe = currentThermalElectron(ne, Te, A)
    return Ithe * (1 + q * (V-Vp) / (k * TeK))**b #nA


def currentIonSaturation(ni, Ti, V, Vp, A, b):
    #[Barjatya +2007 thesis Eqn 1.5]
    #Ion saturation current (NOTE: ion current dominates e- current here)
    A = A * 1e-4 #convert from cm^2 to m^2
    ni = ni * 1e6 #convert from cm^-3 to m^-3
    TiK = Ti * 11604.5 #convert from eV to K
    Ithi = currentThermalIon(ni, Ti, A)
    return Ithi * (1 + q * (V-Vp) / (k * TiK))**b #nA



#------------------------------------------------------
#Drifting, unmagnetized, Maxwellian plasmas, collisionless
#------------------------------------------------------

def currentMesothermalElectronSaturation(ne, Te, V, Vp, A):
    #Mesothermal plasma - relevant to fast-moving satellite where velocity 
    #exceeds ion thermal velocity but not electron thermal velocity
    #Barjatya+2007 thesis Eqn 1.6
    A = A * 1e-4 #convert from cm^2 to m^2
    ne = ne * 1e6 #convert from cm^-3 to m^-3
    TeK = Te * 11604.5 #convert from eV to K
    Ithe = currentThermalElectron(ne, Te, A)
    return Ithe * (me * V**2 / (2 * k * TeK) + 0.5 + q * (V-Vp) / (k * TeK))**0.5 #nA

def currentMesothermalIonSaturation(ni, Ti, V, Vp, A):
    #Barjatya+2007 thesis Eqn 1.6
    A = A * 1e-4 #convert from cm^2 to m^2
    ni = ni * 1e6 #convert from cm^-3 to m^-3
    TiK = Ti * 11604.5 #convert from eV to K
    Ithi = currentThermalIon(ni, Ti, A)
    return Ithi * (mi*V**2 / (2 * k * TiK) + 0.5 + q * (V-Vp) / (k * TiK))**0.5 #nA


#------------------------------------------------------
#Magnetized, non-drifting, nonthermal plasmas, collisionless
#------------------------------------------------------

def currentMagnetizedNonthermalElectronSaturation(ne, Te, V, Vp, B, A):
    #Barjatya+2007 thesis Eqn 1.7 via Parker and Murphy+67: 
    #(NOTE: "This effectively gives a canonical upper bound to the collected saturation current")
    # Parker and Murphy (1967), “Potential buildup on an electron-emitting ionosphericsatelite,” JGR, vol. 72, p. 1631, 1967.
    A = A * 1e-4 #convert from cm^2 to m^2
    ne = ne * 1e6 #convert from cm^-3 to m^-3
    Ithe = currentThermalElectron(ne, Te, A)
    return (Ithe/2) * (1 + (8 * np.abs(q * (V - Vp)) / (me * we**2 * r**2))**0.5) #nA

def currentMagnetizedNonthermalIonSaturation(ni, Ti, V, Vp, B, A):
    #Barjatya+2007 thesis Eqn 1.7 via Parker and Murphy+67: 
    #(NOTE: "This effectively gives a canonical upper bound to the collected saturation current")
    # Parker and Murphy (1967), “Potential buildup on an electron-emitting ionosphericsatelite,” JGR, vol. 72, p. 1631, 1967.
    A = A * 1e-4 #convert from cm^2 to m^2
    ni = ni * 1e6 #convert from cm^-3 to m^-3
    Ithi = currentThermalIon(ni, Ti, A)
    return (Ithi/2) * (1 + (8 * np.abs(q * (V - Vp)) / (mi * wi**2 * r**2))**0.5) #nA


#------------------------------------------------------
#Magnetized, non-drifting, Maxwellian plasmas, collisionless
#------------------------------------------------------

def currentMagnetizedThermalElectronSaturation(ne, Te, V, Vp, B, A):
    #Barjatya+2007 thesis Eqn 1.8 via Rubinstein and Laframboise
    #Rubinstein and Laframboise (1982), “Theory of a spherical probe in a collisionlessmagnetoplasma,” Physics of Fluids, vol. 25, p. 1174
    A = A * 1e-4 #convert from cm^2 to m^2
    ne = ne * 1e6 #convert from cm^-3 to m^-3
    TeK = Te * 11604.5 #convert from eV to K
    Ithe = currentThermalElectron(ne, Te, A)
    return Ithe * (0.5 + 0.5*(8*np.abs(q*(V-Vp))/(me*we**2*r**2))**0.5 + k*TeK/(me*we**2*r**2)) #nA

def currentMagnetizedThermalIonSaturation(ni, Ti, V, Vp, B, A):
    #Barjatya+2007 thesis Eqn 1.8 via Rubinstein and Laframboise
    #Rubinstein and Laframboise (1982), “Theory of a spherical probe in a collisionlessmagnetoplasma,” Physics of Fluids, vol. 25, p. 1174
    A = A * 1e-4 #convert from cm^2 to m^2
    ni = ni * 1e6 #convert from cm^-3 to m^-3
    TiK = Ti * 11604.5 #convert from eV to K
    Ithi = currentThermalIon(ni, Ti, A)
    return Ithi * (0.5 + 0.5*(8*np.abs(q*(V-Vp))/(mi*wi**2*r**2))**0.5 + k*TiK/(mi*wi**2*r**2)) #nA



#------------------------------------------------------
#Determine currents knowing densities and temperatures 
#------------------------------------------------------

ne = 5000 #cm^-3
Te = 0.3 #eV
ni = 5000 #cm^-3
Ti = 0.3 #eV

#See Asenovski+18 (https://doi.org/10.1016/j.jastp.2017.12.013) for a discussion of photoemission current densities
Iph0 = 3 #nA/cm^2 
r = 4 #probe radius (cm)

V = 1 #probe potential (defined relative to sc ground)

#------------------------------------------------------
#Choose probe type
#------------------------------------------------------
#b = 0    #flat probe
#b = 0.5  #cylindrical probe
b = 1    #spherical probe




#Debye length 
lambda_D = debye_length(ne, Te)

#photoelectron current
Iph = currentPhotoelectron(Iph0, Te, V, Vp, r)


#Thermal electron and ion current (Barjatya +2007 thesis Eqn 1.2)
Ithe = currentThermalElectron(ne, Te, A)
Ithi = currentThermalIon(ni, Ti, A)


Ie = currentElectronRetardation(ne, Te, V, Vp, A)


#Electron and ion saturation currents
Isat_e = currentElectronSaturation(ne, Te, V, Vp, A, b)
Isat_i = currentIonSaturation(ni, Ti, V ,Vp, A, b)


Isat_eM = currentMesothermalElectronSaturation(ne, Te, V, Vp, A)
Isat_iM = currentMesothermalIonSaturation(ni, Ti, V, Vp, A)


#-------------------------------
#Magnetized plasmas 
#-------------------------------

B = 100 #nT

we = 28 * B #electron gyrofrequency
wi = 28 * B / 1836 #ion gyrofrequency 

r = 0.01 #probe radius (m)


Isat_e = currentMagnetizedNonthermalElectronSaturation(ne, Te, V, Vp, B, A)
Isat_i = currentMagnetizedNonthermalIonSaturation(ni, Ti, V, Vp, B, A)


Ie = currentMagnetizedThermalElectronSaturation(ne, Te, V, Vp, B, A)
Ii = currentMagnetizedThermalIonSaturation(ni, Ti, V, Vp, B, A)







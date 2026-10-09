"""

Determine currents collected by a Langmuir or electric field probe in various plasma regimes
(unmagnetized, magnetized, thermal, nonthermal, mesothermal) using the equations from Barjatya+2007 thesis.

NOTE: sense of current is defined to be the same as the Langmuir Probe standard. Namely, + current away from probe is positive.

Useful reference is Barjatya+2007 thesis

NOTE: Needs testing against published data!!!!!!!!!!!!!!!!!


NOTE (Error): I need to determine how to extend the ion current into positive potentials.
As it is, the only term is the ion saturation current, which by definition cuts off at Vp=0.
However, there should be an ion thermal term that decays with potential that extends into positive potentials.
This should be the equivalent of an "ion retardation region"


"""

import matplotlib.pyplot as plt
import numpy as np

# ------------------------------------------------------
# Define constants
# ------------------------------------------------------

k = 1.38e-23  # Boltzmann constant (m^2 kg s^-2 K^-1)
q = 1.6e-19  # elementary charge (C)
me = 9.11e-31  # electron mass (kg)
mi = 1.67e-27  # proton mass (kg)
eo = 8.85e-12  # permittivity of free space (s^4 A^2 m^-3 kg^-1)

# ------------------------------------------------------
# Choose probe type
# ------------------------------------------------------
# b = 0    #flat probe
# b = 0.5  #cylindrical probe
# b = 1    #spherical probe


# VERIFIED
def debye_length(ne, Te):
    # ne --> cm-3
    # Te --> eV

    TeK = Te * 11604.5  # convert from eV to K
    ne_goo = ne * 1e6  # convert from cm^-3 to m^-3
    # return np.sqrt(eo * k * TeK / (2 * np.pi * me))
    return np.sqrt(eo * k * TeK / (ne_goo * q**2))  # m


def ramCurrent_ion(ni, A, VelSC):
    # Ion ram current in nA (negative sign because ion current flows TO the probe)
    # Note that the electron thermal velocity generally exceeds the ram velocity, even for LEO satellties at 7.6 km/sec,
    # and so we can generally only consider the ion current.

    # ni --> plasma ion density in cm^-3
    # A --> probe surface area in cm^2
    # VelSC --> spacecraft velocity in m/s

    ni_goo = ni * 1e6  # convert from cm^-3 to m^-3
    A_goo = A / 1e4  # convert from cm^2 to m^2

    return -ni_goo * q * A_goo * VelSC * 1e9
    Iram_uA = -ni_rough * q * AreaProbe * VelSC * 1e6


def ramCurrent_electron(ne, A, VelSC):
    # Electron ram current in nA (positive sign for electrons flowing TO the probe)
    # Note that the electron thermal velocity generally exceeds the ram velocity, even for LEO satellties at 7.6 km/sec,
    # and so we can generally only consider the ion current.

    # ne --> plasma electron density in cm^-3
    # A --> probe surface area in cm^2
    # VelSC --> spacecraft velocity in m/s

    ne_goo = ne * 1e6  # convert from cm^-3 to m^-3
    A_goo = A / 1e4  # convert from cm^2 to m^2

    return ne_goo * q * A_goo * VelSC * 1e9


def potentialFloating_NoPhotoelectrons(Te, amu):
    # Floating potential of a probe in a plasma (no photoemission; assumes Te=Ti)
    # Will be about -1 V in F-region ionosphere assuming O+ and Te=1000 K.
    # Te --> eV

    TeK = Te * 11604.5  # convert from eV to K
    return -1 * (k * (TeK / q)) * np.log(amu * mi / me) ** 0.5  # V


def currentPhotoelectron(Iph0, TePhoto, dV, A):
    # Photoelectron current (e.g. Hinteregger+59) (effectively an ion current to probe)
    # A is effective area of the probe surface (cm^2)
    # Te is the energy of photoelectrons (eV)
    # Iph0 is the photoelectron current density at the probe surface (nA/cm^2). Defined as +
    #   This depends on the material properties of the probe and the solar UV flux,
    #   and can be estimated from models or measurements (e.g. Fahleson+74).

    # NOTE: photoelectrons leaving probe same sense as ions hitting probe, so we define this as a negative current.

    TeK = TePhoto * 11604.5  # convert from eV to K

    if dV < 0:
        return -A * Iph0
    else:
        return -A * Iph0 * np.exp((-q * dV) / (k * TeK))  # nA


# ------------------------------------------------------
# Non-drifting, unmagnetized, Maxwellian plasmas, collisionless
# ------------------------------------------------------


def currentThermalElectron(ne, Te, A):
    # [Barjatya +2007 thesis Eqn 1.2]
    # Electron current to probe defined as +
    # A --> cm^2
    # ne --> cm^-3
    # Te --> eV
    A_goo = A * 1e-4  # convert from cm^2 to m^2
    ne_goo = ne * 1e6  # convert from cm^-3 to m^-3
    TeK = Te * 11604.5  # convert from eV to K
    return 1e9 * ne_goo * q * A_goo * np.sqrt((k * TeK) / (2 * np.pi * me))  # nA


def currentThermalIon(ni, Ti, A, amu):
    # [Barjatya +2007 thesis Eqn 1.2]
    # Ion current to probe defined as -
    # A --> cm^2
    # ni --> cm^-3
    # Ti --> eV
    A_goo = A * 1e-4  # convert from cm^2 to m^2
    ni_goo = ni * 1e6  # convert from cm^-3 to m^-3
    TiK = Ti * 11604.5  # convert from eV to K
    return -1e9 * ni_goo * q * A_goo * np.sqrt((k * TiK) / (2 * np.pi * amu * mi))  # nA


def currentElectrons_ElectronRetardation(ne, Te, dV, A):
    # [Barjatya +2007 thesis Eqn 1.4]
    # Electron current in electron retardation region current (NOTE: e- current dominates ion current here)
    # A --> cm^2
    # ne --> cm^-3
    # Te --> eV
    # dV --> V
    Ithe = currentThermalElectron(ne, Te, A)
    TeK = Te * 11604.5  # convert from eV to K
    if dV < 0:
        return Ithe * np.exp(q * dV / (k * TeK))  # nA
    if dV > 0:
        return np.nan


# ...IS THIS CORRECT?
def currentIons_ElectronRetardation(ni, Ti, dV, A, amu):
    # Ion current in electron retardation region
    # A --> cm^2
    # ni --> cm^-3
    # Ti --> eV
    # dV --> V
    Ithi = currentThermalIon(ni, Ti, A, amu)
    TiK = Ti * 11604.5  # convert from eV to K

    if dV > 0:
        return Ithi * np.exp(-q * dV / (k * TiK))  # nA
    if dV < 0:
        return np.nan


def currentElectronSaturation(ne, Te, dV, A, betaVal):
    # [Barjatya +2007 thesis Eqn 1.5]
    # NOTE: this expression is derived for positive potentials (e.g. see Allen+92 eqns 13, 14)
    # Electron saturation current (NOTE: e- current dominates ion current here)
    # A --> cm^2
    # ne --> cm^-3
    # Te --> eV
    # dV --> V
    Ithe = currentThermalElectron(ne, Te, A)
    TeK = Te * 11604.5  # convert from eV to K

    if dV < 0:
        return np.nan
    if dV > 0:
        return Ithe * (1 + (q * dV / (k * TeK))) ** betaVal  # nA


def currentIonSaturation(ni, Ti, dV, A, betaVal, amu):
    # [Barjatya +2007 thesis Eqn 1.5]
    # NOTE: this expression is derived for negative potentials (e.g. see Allen+92 eqn 17)
    # Ion saturation current (NOTE: ion current dominates e- current here)
    # A --> cm^2
    # ni --> cm^-3
    # Ti --> eV
    # dV --> V
    Ithi = currentThermalIon(ni, Ti, A, amu)
    TiK = Ti * 11604.5  # convert from eV to K

    if dV < 0:
        return Ithi * (1 + (-q) * dV / (k * TiK)) ** betaVal  # nA
    if dV > 0:
        return np.nan


# ------------------------------------------------------
# Drifting, unmagnetized, Maxwellian plasmas, collisionless
# ------------------------------------------------------


def currentMesothermalElectronSaturation(ne, Te, vel, dV, A):
    # Mesothermal plasma - relevant to fast-moving satellite where velocity
    # exceeds ion thermal velocity but not electron thermal velocity
    # Barjatya+2007 thesis Eqn 1.6
    # A --> cm^2
    # ne --> cm^-3
    # Te --> eV
    # dV --> V
    # vel --> m/s
    Ithe = currentThermalElectron(ne, Te, A)
    TeK = Te * 11604.5  # convert from eV to K
    t1 = me * vel**2 / (2 * k * TeK)
    t2 = 0.5
    t3 = q * dV / (k * TeK)
    # t3 = q * (V - Vp) / (k * TeK)
    return Ithe * (t1 + t2 + t3) ** 0.5  # nA


def currentMesothermalIonSaturation(ni, Ti, vel, dV, A, amu):
    # Barjatya+2007 thesis Eqn 1.6
    # A --> cm^2
    # ni --> cm^-3
    # Ti --> eV
    # dV --> V
    # vel --> m/s
    Ithi = currentThermalIon(ni, Ti, A, amu)
    TiK = Ti * 11604.5  # convert from eV to K
    t1 = amu * mi * vel**2 / (2 * k * TiK)
    t2 = 0.5
    t3 = -q * dV / (k * TiK)
    # t3 = q * (V - Vp) / (k * TiK)
    return Ithi * (t1 + t2 + t3) ** 0.5  # nA


# ------------------------------------------------------
# Magnetized, non-drifting, nonthermal plasmas, collisionless
# ------------------------------------------------------


def currentMagnetizedNonthermalElectronSaturation(ne, Te, dV, B, A, r):
    # Barjatya+2007 thesis Eqn 1.7 via Parker and Murphy+67:
    # NOTE: spherical probe only
    # (NOTE: "This effectively gives a canonical upper bound to the collected saturation current")
    # Parker and Murphy (1967), “Potential buildup on an electron-emitting ionosphericsatelite,” JGR, vol. 72, p. 1631, 1967.
    #   B = magnetic field strength (nT)
    #   r = probe radius (cm)
    #   Te = electron temperature (eV)
    #   A = probe surface area (cm^2)
    #   dV --> V
    #   ne --> cm^-3
    Ithe = currentThermalElectron(ne, Te, A)
    wc = 2 * np.pi * 28 * B  # electron gyrofrequency (rad/s)
    rm = r * 0.01  # m
    return (Ithe / 2) * (1 + (8 * np.abs(q * dV) / (me * wc**2 * rm**2)) ** 0.5)  # nA


def currentMagnetizedNonthermalIonSaturation(ni, Ti, dV, B, A, r, amu):
    # Barjatya+2007 thesis Eqn 1.7 via Parker and Murphy+67:
    # NOTE: spherical probe only
    # (NOTE: "This effectively gives a canonical upper bound to the collected saturation current")
    # Parker and Murphy (1967), “Potential buildup on an electron-emitting ionosphericsatelite,” JGR, vol. 72, p. 1631, 1967.
    # A --> cm^2
    # ni --> cm^-3
    # Ti --> eV
    # dV --> V
    # B --> nT
    # r --> cm
    Ithi = currentThermalIon(ni, Ti, A, amu)
    wc = 2 * np.pi * 28 * B / (amu * 1836)
    rm = r * 0.01  # m
    return (Ithi / 2) * (
        1 + (8 * np.abs(q * dV) / (amu * mi * wc**2 * rm**2)) ** 0.5
    )  # nA


# ------------------------------------------------------
# Magnetized, non-drifting, Maxwellian plasmas, collisionless
# ------------------------------------------------------


def currentMagnetizedThermalElectronSaturation(ne, Te, dV, B, A, r):
    # Barjatya+2007 thesis Eqn 1.8 via Rubinstein and Laframboise
    # NOTE: spherical probe only
    # Rubinstein and Laframboise (1982), “Theory of a spherical probe in a collisionlessmagnetoplasma,” Physics of Fluids, vol. 25, p. 1174
    # A --> cm^2
    # ne --> cm^-3
    # Te --> eV
    # dV --> V
    # B --> nT
    # r --> cm
    Ithe = currentThermalElectron(ne, Te, A)
    TeK = Te * 11604.5  # convert from eV to K
    wc = 2 * np.pi * 28 * B
    rm = r * 0.01  # m

    t1 = 0.5
    t2 = 0.5 * (8 * np.abs(q * dV) / (me * wc**2 * rm**2)) ** 0.5
    # t2 = 0.5 * (8 * np.abs(q * (V - Vp)) / (me * wc**2 * rm**2))**0.5
    t3 = k * TeK / (me * wc**2 * rm**2)

    return Ithe * (t1 + t2 + t3)  # nA


def currentMagnetizedThermalIonSaturation(ni, Ti, dV, B, A, r, amu):
    # Barjatya+2007 thesis Eqn 1.8 via Rubinstein and Laframboise
    # NOTE: spherical probe only
    # Rubinstein and Laframboise (1982), “Theory of a spherical probe in a collisionlessmagnetoplasma,” Physics of Fluids, vol. 25, p. 1174
    # A --> cm^2
    # ni --> cm^-3
    # Ti --> eV
    # dV --> V
    # B --> nT
    # r --> cm
    Ithi = currentThermalIon(ni, Ti, A, amu)
    TiK = Ti * 11604.5  # convert from eV to K
    wc = 2 * np.pi * 28 * B / (amu * 1836)  # electron gyrofrequency (rad/sec)
    rm = r * 0.01  # m
    t1 = 0.5
    t2 = 0.5 * (8 * np.abs(q * dV) / (amu * mi * wc**2 * rm**2)) ** 0.5
    # t2 = 0.5 * (8 * np.abs(q * (V - Vp)) / (amu * mi * wc**2 * rm**2))**0.5
    t3 = k * TiK / (amu * mi * wc**2 * rm**2)

    return Ithi * (t1 + t2 + t3)  # nA

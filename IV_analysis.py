# **********************************
"""
Simulate expected I-V curves based on instrument Langmuir probe OR Electric field probe specs

*Use to estimate max and min currents for amplifier selection
*Use to estimate bias currents required to keep the probe in the correct operating regime
*Use to estimate sampling resolution required for ion and electron saturation regions, and electron retardation region

NOTE: this is the generalized routine. Make copies of this and adapt it to specific missions.

NOTE: I'm not considering secondary electrons here. These can be very important for dealing with high energy electron flux.
See Ch. 5 of Lai and Pradipta book Physics of Satellite Surface Charging.


"""

import sys

sys.path.append("/Users/abrenema/Desktop/code/Aaron/github/antenna_analysis/")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from langmuirProbes import *

# **********************************

payload = "TBD Mission"  # name of mission

eclipse = False

inststr = "Cylindrical"
# inststr = "Planar"


if inststr == "Cylindrical":
    # SLP Probe properties
    betaVal = 0.5  # cylindrical probe
    r = 0.32  # cm (cylindrical probe radius)
    l = 5  # cm (probe length)
    ## Effective area of cylindrical probe - i.e. rectangular cross section.
    Aeff = 2 * r * l


# Fixed bias probe (patch probe)
if inststr == "Planar":
    betaVal = 0  # planar probe
    r = 3.95 / 2.0  # cm (planar probe radius)
    ## r = 6 / 2.0  # cm (planar probe radius)
    Aeff = np.pi * r**2  # cm^2 (probe area)


# C/NOFS example: 100 pA for a few*10^2/cc for a 2.6 cm diameter sphere and -5V bias
# betaVal = 1  # sphere
# r = 2.6 / 2  # cm (spherical probe radius)
# Aeff = np.pi * r**2  # sunlit area only


B = 50000  # nT
VelSC = 7600  # m/sec


# Range of densities we're required to work over is 10^2 - 10^7 cm-3
ne = 1e1  # cm^-3
ni = ne

# Range of electron temps from 400-10000 K (0.03-0.86 eV)
# Te = 0.86  # eV
Te = 4  # eV
# Te = 0.03  # eV
# Te = 0.3
Ti = Te
# Ti = 0.05  # eV
amu = 1

# See Asenovski+18 (https://doi.org/10.1016/j.jastp.2017.12.013) for a discussion of photoemission current densities
Iph0 = 3  # nA/cm^2 (current density estimate for Gold [e.g. Hartzell+23]). Appropriate for ionosphere
# Iph0 = 6  # nA/cm^2 (estimate for TiN [Laakso+26; Pederson+08]). Appropriate for solar wind and low density magnetosphere

# Photoelectron temperature (eV)
TePhoto = 5  # (1-10 eV typically)


# add in a bias current in nA (e.g. 70% of the max photoelectron current)
# NOTE: that we usually bias negative currents. However, here an electron current to the probe is defined as positive.
Ibias = 0.7 * 0.009 * 1000
# Ibias = 0

# Debye length (few cm)
lambda_D = debye_length(ne, Te)  # m


Vf = potentialFloating_NoPhotoelectrons(Te, amu)


# Thermal electron and ion current (Barjatya +2007 thesis Eqn 1.2)
Ithe = currentThermalElectron(ne, Te, Aeff)
Ithi = currentThermalIon(ni, Ti, Aeff, amu)


# dVr = np.arange(-6, 15, 0.01)  # probe potential relative to plasma potential (V - Vp)
dVr = np.arange(-6, 15, 0.01)  # probe potential relative to plasma potential (V - Vp)

Iph = np.empty(len(dVr))
Ie_ret = np.empty(len(dVr))
Ii_ret = np.empty(len(dVr))
Ie_sat = np.empty(len(dVr))
Ii_sat = np.empty(len(dVr))
Ie_satMeso = np.empty(len(dVr))
Ii_satMeso = np.empty(len(dVr))
Ii_ram = np.empty(len(dVr))
Ibias = np.full(len(dVr), Ibias)


for i in range(len(dVr)):
    # Relevant b/t Vf and Vp
    Ie_ret[i] = currentElectrons_ElectronRetardation(ne, Te, dVr[i], Aeff)
    Ii_ret[i] = currentIons_ElectronRetardation(ni, Ti, dVr[i], Aeff, amu)

    # Photoelectron current
    Iph[i] = currentPhotoelectron(Iph0, TePhoto, dVr[i], Aeff)  # nA

    # Electron and ion saturation currents
    Ie_sat[i] = currentElectronSaturation(ne, Te, dVr[i], Aeff, betaVal)
    Ie_satMeso[i] = currentMesothermalElectronSaturation(ne, Te, VelSC, dVr[i], Aeff)

    Ii_sat[i] = currentIonSaturation(ni, Ti, dVr[i], Aeff, betaVal, amu)
    Ii_satMeso[i] = currentMesothermalIonSaturation(ni, Ti, VelSC, dVr[i], Aeff, amu)

    Ii_ram[i] = ramCurrent_ion(ni, Aeff, VelSC)  # nA


# -----------------------------------------------------------------------------------------------------

if eclipse:
    print("Eclipse: photoemission current is zero")
    Iph = np.zeros(len(dVr))


# sum all the currents
Itot = np.empty(len(dVr))
for i in range(len(Itot)):
    Itot[i] = np.nansum(
        [Ie_ret[i], Ie_sat[i], Ii_sat[i], Ii_ret[i], Iph[i], Ii_ram[i], Ibias[i]]
    )


# Change color of total current based on which sub current dominates
current_magnitudes = np.column_stack(
    [
        np.abs(Ie_ret),
        np.abs(Ie_sat),
        np.abs(Ii_sat),
        np.abs(Ii_ret),
        np.abs(Iph),
        np.abs(Ii_ram),
        np.abs(Ibias),
    ]
)
dominant_current = np.nanargmax(current_magnitudes, axis=1)
dominant_colors = np.array(
    ["green", "red", "blue", "blue", "goldenrod", "purple", "darkturquoise"]
)[dominant_current]


itot_values = np.abs(Itot) / 1000
itot_points = np.column_stack([dVr, itot_values])
itot_segments = np.stack([itot_points[:-1], itot_points[1:]], axis=1)
itot_line = LineCollection(itot_segments, colors=dominant_colors[:-1], linewidths=1.5)


# Identify the value of the floating potential (relative to the plasma potential at zero)
Vf = dVr[np.nanargmin(itot_values)]


# Find max/min values for some of the plots.
minv = np.nanmin([Ie_ret, Ie_sat, Ii_ram, Iph, Ii_sat, Itot, Ibias]) / 1000
maxv = np.nanmax([Ie_ret, Ie_sat, Ii_ram, Iph, Ii_sat, Itot, Ii_ret, Ibias]) / 1000
minvp = np.nanmin([Ie_ret, Ie_sat]) / 1000
maxvp = np.nanmax([Ie_ret, Ie_sat]) / 1000
minvm = np.nanmin([Ii_ram, Iph, Ii_sat, Ii_ret]) / 1000
maxvm = np.nanmax([Ii_ram, Iph, Ii_sat, Ii_ret]) / 1000


# Plot individual currents
fig, axs = plt.subplots(4, figsize=(9, 9))
fig.suptitle(
    "I-V curves for "
    + payload
    + " "
    + inststr
    + " probe"
    + "\nBo="
    + str(B)
    + " nT"
    + "; VelSC="
    + str(VelSC)
    + " m/s"
    + "\nne="
    + f"{ne:.2e}"
    + " cm-3; ni="
    + f"{ni:.2e}"
    + " cm-3\nTe="
    + str(Te)
    + " eV; Ti="
    + f"{Ti:.2g}"
    + " eV; amu="
    + str(amu)
    + "\neclipse="
    + str(eclipse)
)
axs[0].plot(dVr, Ie_ret / 1000, label="Ie_ret", color="green")
axs[0].set_ylabel("Ie_ret (uA)")
axs[0].plot(dVr, Ie_sat / 1000, label="Ie_sat", color="red")
axs[0].set_ylabel("+ currents (uA)")
if np.nansum(Ibias) > 0:
    axs[0].plot(dVr, Ibias / 1000, label="Ibias", color="darkturquoise", linestyle="-.")
axs[1].plot(dVr, np.abs(Ii_sat) / 1000, label="|Ii_sat|", color="blue")
axs[1].plot(dVr, np.abs(Ii_ret) / 1000, label="|Ii_ret|", color="blue")
axs[1].set_ylabel("-Ii_sat (uA)")
axs[1].plot(dVr, np.abs(Iph) / 1000, label="|Iph|", color="goldenrod")
axs[1].set_ylabel("-Iph (uA)")
axs[1].plot(dVr, np.abs(Ii_ram) / 1000, label="|Ii_ram|", color="purple")
axs[1].set_ylabel("- currents (uA)")
if np.nansum(Ibias) < 0:
    axs[1].plot(
        dVr, np.abs(Ibias) / 1000, label="Ibias", color="darkturquoise", linestyle="-."
    )
axs[2].add_collection(itot_line)
axs[2].plot([], [], color="green", label="Ie_ret")
axs[2].plot([], [], color="red", label="Ie_sat")
axs[2].plot([], [], color="blue", label="|Ii_sat|")
axs[2].plot([], [], color="blue", label="|Ii_ret|")
axs[2].plot([], [], color="goldenrod", label="|Iph|")
axs[2].plot([], [], color="purple", label="|Ii_ram|")
axs[2].plot([], [], color="darkturquoise", label="|Ibias|")
axs[2].set_ylabel("Itot (uA)")
axs[2].set_xlabel("volts")
axs[3].plot(dVr, Ie_ret / 1000, color="green", label="Ie_ret")
axs[3].plot(dVr, Ie_sat / 1000, color="red", label="Ie_sat")
axs[3].plot(dVr, Ii_ram / 1000, color="purple", label="Ii_ram")
axs[3].plot(dVr, Iph / 1000, color="goldenrod", label="Iph")
axs[3].plot(dVr, Ii_sat / 1000, color="blue", label="Ii_sat")
axs[3].plot(dVr, Ii_ret / 1000, color="blue", label="Ii_ret")
axs[3].plot(dVr, Itot / 1000, color="black", linestyle="--", linewidth=2, label="Itot")
axs[3].plot(
    dVr, Ibias / 1000, color="darkturquoise", linestyle="-.", linewidth=2, label="Ibias"
)
for i in range(4):
    axs[i].set_xlim(np.min(dVr), np.max(dVr))
    axs[i].set_yscale("log")
    axs[i].legend()
    axs[i].axvline(
        x=Vf, color="red", linestyle="--", linewidth=0.5, label="Floating Potential"
    )
    axs[i].axvline(
        x=0, color="blue", linestyle="--", linewidth=0.5, label="Plasma Potential"
    )
axs[0].set_ylim(minvp / 100, maxvp * 100)
axs[1].set_ylim(minvm / 100, maxvm * 100)
axs[2].set_ylim(minv / 100, maxv * 100)
axs[3].set_yscale("linear")
axs[3].set_ylim(minv + 0.4 * minv, maxv + 0.4 * maxv)
axs[3].set_ylabel("Currents (uA)\nLinear scale")
axs[3].set_xlabel("Probe Potential Relative to Plasma Potential")


#############################################################################
# The below code is not currently working


# Save total current as pickle file so I can use it in the TBLZ_langmuirProbe_IV_extract.py script to test these curves.
# i.e. make sure that I'm extracting the same ne, Te, ni, Vp that I'm inputting.
import pickle

with open(
    "/Users/abrenema/Desktop/code/Aaron/github/mission_routines/trailblazer/TBLZ_LangmuirProbe_Itot.pkl",
    "wb",
) as f:
    pickle.dump((Itot, dVr), f)


# tst = np.where(dVr > 0)[0][0]


# Relevant b/t Vf and Vp
Ie_ret = currentElectrons_ElectronRetardation(ne, Te, dV, Aeff)


# Electron and ion saturation currents
Ie_sat = currentElectronSaturation(ne, Te, dV, Aeff, betaVal)
Ie_satMeso = currentMesothermalElectronSaturation(ne, Te, VelSC, dV, Aeff)

Ii_sat = currentIonSaturation(ni, Ti, dV, Aeff, betaVal, amu)
Ii_satMeso = currentMesothermalIonSaturation(ni, Ti, VelSC, dV, Aeff, amu)

# --------------------------------------------------------------------------------
# Find voltage limits required for a range of temperatures and densities.


dVr = np.arange(-1, 1, 0.05)


Ie_ret400k = np.empty(len(dVr))
Ie_ret10000k = np.empty(len(dVr))
for i in range(len(dVr)):
    Ie_ret400k[i] = currentElectrons_ElectronRetardation(4000, 0.03, dVr[i], Aeff)
    Ie_ret10000k[i] = currentElectrons_ElectronRetardation(4000, 0.86, dVr[i], Aeff)


# Plot individual currents
fig, axs = plt.subplots(2)
axs[0].plot(dVr, Ie_ret400k)
axs[0].set_ylim(0, np.nanmax(Ie_ret400k))
axs[1].plot(dVr, Ie_ret10000k)
axs[1].set_ylim(0, np.nanmax(Ie_ret10000k))

for i in range(2):
    axs[i].set_xlim(-1.5, 1.5)
    axs[i].set_yscale("log")
    axs[i].set_ylim(0.1, 1e14)

    # axs[i].set_xlim(np.nanmin(dVr), np.nanmax(dVr))

plt.show()
plt.show()

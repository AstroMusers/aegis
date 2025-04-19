import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import pandas as pd
import sys

from radio_module import *


# Matplotlib Configuration
# ------------------------
mpl.rcParams["figure.autolayout"] = True
plt.rcParams['figure.figsize'] = [8, 4]

rc = {"font.family": "times new roman",
      "font.size": 11,
      "mathtext.fontset": "stix"}
plt.rcParams.update(rc)

mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=["firebrick", "navy", "b"])

# Random number generator
rng = np.random.default_rng()

# Get the input variables from command-line arguments. Use this if you want
# to run the code using the GUI, otherwise comment out the following
# 11 lines.

a_min = float(sys.argv[1])
a_max = float(sys.argv[2])
B_mean = float(sys.argv[3])
B_sd = float(sys.argv[4])
M_mean = float(sys.argv[5])
M_sd = float(sys.argv[6])
Mdot_mean = float(sys.argv[7])
Mdot_sd = float(sys.argv[8])
D_mean = float(sys.argv[9])
D_sd = float(sys.argv[10])
size = int(sys.argv[11])
IsBurst = int(sys.argv[12])

# Default values for the variables. Use this section of the code to run
# the code in the script.
# --------------------------------------------------------------------
# size = 150
# a_min = -2
# a_max = 2
# B_mean = 10
# B_sd = 3
# M_mean = 3
# M_sd = 0.8
# Mdot_mean = 100
# Mdot_sd = 30
# D_mean = 1000
# D_sd = 300
# IsBurst = True

# --------------------------------------------------------------------

# Create Synthetic Sample
# -----------------------
exoplanets = []
for i in range(size):
    a = rng.uniform(a_min, a_max)  # logAU

    Bs = rng.normal(B_mean, B_sd)  # gauss
    assert Bs >= 0, f"Magnetic Field Strength is expected nonnegative, instead got {Bs=}."

    M = rng.normal(M_mean, M_sd)  # solar masses
    assert M > 0, f"Stellar Mass is expected positive, instead got {M=}."

    Mdot = rng.normal(Mdot_mean, Mdot_sd)  # e-15 solar masses / yr
    assert Mdot > 0, f"Stellar Mass Loss Rate is expected positive, instead got {Mdot=}."

    D = rng.normal(D_mean, D_sd)  # light years
    assert D > 0, f"Distance to Exoplanet is expected positive, instead got {D=}."

    exo = Exoplanet(i, a, 1, 1, 1, Bs, M, Mdot, D)  # See the module. 1, 1, 1 stand for currently unnecessary parameters.
    exoplanets.append(exo)

print(exoplanets)

# Currently unused parameters
v_sw = 400  # km/s
B_sw = 10 ** (-8)  # teslas

# Make Calculations
# -----------------
frequencies = []
intensities = []
distances = []
semis = []

for exo in exoplanets:
    a = exo.semi_major_axis
    B = exo.magnetic_field
    M = exo.star_mass
    Mdot = exo.star_mass_loss
    D = exo.distance

    semis.append(a)
    distances.append(D)

    a = 10 ** a  # conversion to AU
    D = D * 9.46 * 10 ** 15  # conversion to meters

    nu = max_freq(B)
    assert nu > 0, f"Maximum emission frequency must be positive, instead got {nu=}."
    frequencies.append(nu/10**6)

    I = complete(B, a, M, Mdot, D)
    assert I > 0, f"Radio brightness must be positive, instead got {I=}."
    intensities.append(I)

frequencies = np.array(frequencies)

distances = np.array(distances)
distances = np.reciprocal(distances)
distances *= 10 ** 4

intensities = np.array(intensities)
burst = intensities * (10 ** 1.53)

semis = np.array(semis)

# Import Ashtari2022 results
df2 = pd.read_csv("dataAshtari22.txt", delimiter="\t", names=["Planet", "freq", "lowSflux", "highSflux"])
df2.freq /= 10**6
df2.lowSflux *= 10**1.53
df2.highSflux *= 10**1.53

if IsBurst:
    df = pd.DataFrame({"x": frequencies,
                       "y": burst,
                       "d": distances,
                       "s": semis})
else:
    df = pd.DataFrame({"x": frequencies,
                       "y": intensities,
                       "d": distances,
                       "s": semis})

# Plotting
# -------------------

fig, ax = plt.subplots()

# Uncomment these to see the random sample predictions
im = ax.scatter(df.x, df.y, c=df.s, s=df.d, cmap="jet_r")
fig.colorbar(im, ax=ax, label="Distance to Host Star ($\log_{10}{\mathrm{(AU)}}$)")

# im = ax.scatter(df2.freq, df2.highSflux, s=10, marker="v")

plt.axvline(x=10, color="black", linestyle="dashed")

ax.set_xscale("log")
ax.set_yscale("log")
# ax.set_xlim(left=10**(-3))
# ax.set_xlim(10**(-2), 10**5)


ax.axvspan(0, 10, alpha=0.2, color="teal")

if IsBurst:
    ax.set_title("Frequency and Intensity of Burst CMI Emissions of a Synthetic Exoplanet Sample")
else:
    ax.set_title("Frequency and Intensity of Quiescent CMI Emissions of a Synthetic Exoplanet Sample")
ax.set_xlabel("Emission Frequency (MHz)")
ax.set_ylabel("Radio Brightness (Jy)")
plt.show()

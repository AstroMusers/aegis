import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import pandas as pd

from radio_module import *

rng = np.random.default_rng()

mpl.rcParams["figure.autolayout"] = True
plt.rcParams['figure.figsize'] = [8, 4]

rc = {"font.family": "times new roman",
      "font.size": 11,
      "mathtext.fontset": "stix"}
plt.rcParams.update(rc)

df = pd.read_csv("NASA2207.csv", skiprows=34)
print(df)

N = 100

totalFreq = np.zeros(len(df.index))
totalI = np.zeros(len(df.index))

for j in range(N):

    exoplanets = []

    for i in df.iterrows():
        name = i[1][0]
        a = i[1][3]
        M = i[1][7]
        t = i[1][11]
        d = i[1][15]

        d *= 3.26156

        highS_Mdot = t ** (-1.23) * 10 ** 3
        lowS_Mdot = t ** (-0.9) * 10 ** 3

        B = rng.uniform(5,15)

        exo = Exoplanet(name, a, B, M, highS_Mdot, d)
        exoplanets.append(exo)

    # print(exoplanets)

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

        a = np.log10(a)
        semis.append(a)
        a = 10 ** a

        distances.append(D)

        D = D * 9.46 * 10 ** 15  # conversion to meters

        nu = max_freq(B)
        assert nu > 0, f"Maximum emission frequency must be positive, instead got {nu=}."
        frequencies.append(nu / 10 ** 6)

        I = complete(B, a, M, Mdot, D)
        assert I > 0, f"Radio brightness must be positive, instead got {I=}."
        intensities.append(I)

    frequencies = np.array(frequencies)
    intensities = np.array(intensities)

    totalFreq += frequencies
    totalI += intensities

frequencies = totalFreq / N
intensities = totalI / N

distances = np.array(distances)
distances = np.reciprocal(distances)
distances *= 10 ** 2.5

intensities = np.array(intensities)
burst = intensities * (10 ** 1.53)

semis = np.array(semis)

IsBurst = 1

if IsBurst:
    df1 = pd.DataFrame({"x": frequencies,
                        "y": burst,
                        "d": distances,
                        "s": semis})
else:
    df1 = pd.DataFrame({"x": frequencies,
                        "y": intensities,
                        "d": distances,
                        "s": semis})

fig0, ax0 = plt.subplots()

# Uncomment these to see the random sample predictions
im = ax0.scatter(df1.x, df1.y, c=df1.s, s=df1.d, cmap="jet_r")
fig0.colorbar(im, ax=ax0, label="Distance to Host Star ($\log_{10}{\mathrm{(AU)}}$)")

# im = ax.scatter(df2.freq, df2.highSflux, s=10, marker="v")

ax0.axvline(x=10, color="black", linestyle="dashed")

ax0.set_xscale("log")
ax0.set_yscale("log")
ax0.set_xlim(6, 30)
# ax.set_xlim(10**(-2), 10**5)


ax0.axvspan(0, 10, alpha=0.2, color="teal")
fig0.tight_layout()

if IsBurst:
    ax0.set_title("Frequency and Intensity of Burst CMI Emissions of the Exoplanet Sample")
else:
    ax0.set_title("Frequency and Intensity of Quiescent CMI Emissions of the Exoplanet Sample")
ax0.set_xlabel("Emission Frequency (MHz)")
ax0.set_ylabel("Radio Brightness (Jy)")

TheBins = np.logspace(-6, 3, 10)

plt.rcParams['figure.figsize'] = [6, 4]

fig1, ax1 = plt.subplots()

if IsBurst:
    ax1.hist(burst, bins=TheBins, edgecolor="black")
    ax1.set_xlabel("Intensity of Burst Emission (Jy)")
    ax1.set_title("Histogram of Burst Emission Intensities")
else:
    ax1.hist(intensities, bins=TheBins, edgecolor="black")
    ax1.set_xlabel("Intensity of Quiescent Emission (Jy)")
    ax1.set_title("Histogram of Quiescent Emission Intensities")

ax1.set_xscale("log")
ax1.set_ylabel("Number of Exoplanets")

plt.show()

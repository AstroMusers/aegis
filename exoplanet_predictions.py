import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import pandas as pd

from radio_module import *

df = pd.read_csv("NASA2207.csv", skiprows=34)
print(df)

exoplanets = []

for i in df.iterrows():
    name = i[1][0]
    a = i[1][3]
    M = i[1][7]
    t = i[1][11]
    d = i[1][15]

    d *= 3.26156

    highS_Mdot = t ** (-1.23) * 10**3
    lowS_Mdot = t ** (-0.9) * 10**3

    exo = Exoplanet(name, a, 1, M, highS_Mdot, d)
    exoplanets.append(exo)

print(exoplanets)

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

distances = np.array(distances)
distances = np.reciprocal(distances)
distances *= 10 ** 2

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

fig, ax = plt.subplots()

# Uncomment these to see the random sample predictions
im = ax.scatter(df1.x, df1.y, c=df1.s, s=df1.d, cmap="jet_r")
fig.colorbar(im, ax=ax, label="Distance to Host Star ($\log_{10}{\mathrm{(AU)}}$)")

# im = ax.scatter(df2.freq, df2.highSflux, s=10, marker="v")

plt.axvline(x=10, color="black", linestyle="dashed")

ax.set_xscale("log")
ax.set_yscale("log")
# ax.set_xlim(6, 30)
# ax.set_xlim(10**(-2), 10**5)


ax.axvspan(0, 10, alpha=0.2, color="teal")

if IsBurst:
    ax.set_title("Frequency and Intensity of Burst CMI Emissions of a Synthetic Exoplanet Sample")
else:
    ax.set_title("Frequency and Intensity of Quiescent CMI Emissions of a Synthetic Exoplanet Sample")
ax.set_xlabel("Emission Frequency (MHz)")
ax.set_ylabel("Radio Brightness (Jy)")
plt.show()

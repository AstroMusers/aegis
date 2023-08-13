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

df = pd.read_csv("NASA0808.csv", skiprows=55)
print(df)

exoplanets = []

for i in df.iterrows():
    j = i[1]
    name = j[0]
    T = j[2]
    a = j[6]
    R = j[10]
    if j[14] == "Mass":
        M = j[14]
    else:
        M = j[14] * 1.15  # Expected Value of the mass based on projected mass
    p = j[19]
    M_s = j[28]
    t = j[32]
    d = j[36]

    p_c = density(p)
    w_p = rotation(T, a)

    d *= 3.26156

    highS_Mdot = t ** (-1.23) * 10 ** 3
    lowS_Mdot = t ** (-0.9) * 10 ** 3

    Mdot = 8.1 * t**(-1.37)

    B = 1  # Tentative
    sigma = 1  # Jupiter conductivity

    exo = Exoplanet(name, a, R, M, p, B, M_s, Mdot, d)

    r_c = convective_radius(exo)
    mu = magnetic_moment(p_c, w_p, r_c, sigma)

    B = magnetic_field(mu, exo.radius)

    exo.magnetic_field = B

    if exo.magnetic_field != 0:
        exoplanets.append(exo)

print(exoplanets)

frequencies = []
intensities = []
distances = []
semis = []
magnetic_fields = []

for exo in exoplanets:
    a = exo.semi_major_axis
    B = exo.magnetic_field
    M_s = exo.star_mass
    Mdot = exo.star_mass_loss
    D = exo.distance

    a = np.log10(a)
    semis.append(a)
    a = 10 ** a

    distances.append(D)

    magnetic_fields.append(B)

    D = D * 9.46 * 10 ** 15  # conversion to meters

    nu = max_freq(B)
    assert nu > 0, f"Maximum emission frequency must be positive, instead got {nu=}."
    frequencies.append(nu / 10 ** 6)

    I = complete(B, a, M_s, Mdot, D)
    assert I > 0, f"Radio brightness must be positive, instead got {I=}."
    intensities.append(I)

frequencies = np.array(frequencies)

distances = np.array(distances)
distances = np.reciprocal(distances)
distances *= 10 ** 2.7

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

im = ax0.scatter(df1.x, df1.y, c=df1.s, s=df1.d, cmap="jet_r")
fig0.colorbar(im, ax=ax0, label="Distance to Host Star ($\log_{10}{\mathrm{(AU)}}$)")

ax0.axvline(x=10, color="black", linestyle="dashed")

ax0.set_xscale("log")
ax0.set_yscale("log")
# ax0.set_xlim(6, 30)
ax0.set_xlim(0.5, 400)


ax0.axvspan(0, 10, alpha=0.2, color="teal")
fig0.tight_layout()

if IsBurst:
    ax0.set_title("Frequency and Intensity of Burst CMI Emissions of the Exoplanet Sample")
else:
    ax0.set_title("Frequency and Intensity of Quiescent CMI Emissions of the Exoplanet Sample")
ax0.set_xlabel("Emission Frequency (MHz)")
ax0.set_ylabel("Radio Brightness (Jy)")

TheBins1 = np.logspace(-9, 2, 12)

plt.rcParams['figure.figsize'] = [6, 4]
rc = {"font.family": "times new roman",
      "font.size": 14,
      "mathtext.fontset": "stix"}
plt.rcParams.update(rc)

fig1, ax1 = plt.subplots()

if IsBurst:
    ax1.hist(burst, bins=TheBins1, edgecolor="black")
    ax1.set_xlabel("Intensity of Burst Emission (Jy)")
    ax1.set_title("Histogram of Burst Emission Intensities")
else:
    ax1.hist(intensities, bins=TheBins1, edgecolor="black")
    ax1.set_xlabel("Intensity of Quiescent Emission (Jy)")
    ax1.set_title("Histogram of Quiescent Emission Intensities")

ax1.set_xscale("log")
ax1.set_ylabel("Number of Exoplanets")

fig2, ax2 = plt.subplots()

TheBins2 = np.logspace(-4, 4, 17)

ax2.hist(magnetic_fields, bins=TheBins2, edgecolor="black")
ax2.set_xlabel("Magnetic Field Strength at the Surface (Gauss)")
ax2.set_title("Histogram of the Magnetic Field Strengths")

ax2.set_xscale("log")
ax2.set_ylabel("Number of Exoplanets")

plt.show()

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

lofar = pd.read_csv("sensitivities.csv")  # Obtained from van Haarlem et al. (2017), 8h integration time, 4.66MHz effective bandwidth
lofar.columns = lofar.iloc[0]
lofar = lofar[1:]
lofar = lofar.drop([1, 2, 3, 11, 12, 13])
lofar = lofar.reset_index(drop=True)
lofar = lofar.apply(pd.to_numeric, errors="ignore")
lofar["NL Core"] = lofar["NL Core"].multiply(10**(-3)) * 5  # 5 sigma sensitivity in Jy
lofar["Full EU"] = lofar["Full EU"].multiply(10**(-3)) * 5  # 5 sigma sensitivity in Jy

L_NL = lofar["NL Core"]
L_EU = lofar["Full EU"]
Freq = lofar["Freq."]

d = {"Bands": ["Band 1", "Band 2", "Band 3", "Band 4"],
     "Frequencies": [[120, 250], [250, 500], [550, 850], [1050, 1450]],  # MHz
     "RMS Noise": [np.array([190, 190]), np.array([50, 50]), np.array([40, 40]), np.array([45, 45])]  # microJy, 10min integration time, 100MHz Bandwidth
     }

uGMRT = pd.DataFrame(data=d)
integration_time = 8 * 60  # minutes
bandwidth = 100  # MHZ
uGMRT["RMS Noise"] = uGMRT["RMS Noise"] * (np.sqrt((100 * 10) / (bandwidth * integration_time))) * 10**(-6) * 5  # 5 sigma sensitivity in Jy

df = pd.read_csv("NASA0808.csv", skiprows=55)
# print(df)

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

# print(exoplanets)

frequencies = []
intensities = []
distances = []
semis = []
magnetic_fields = []
labels = []

IsBurst = 1

for exo in exoplanets:
    a = exo.semi_major_axis
    B = exo.magnetic_field
    M_s = exo.star_mass
    Mdot = exo.star_mass_loss
    D = exo.distance
    obs = False

    a = np.log10(a)
    semis.append(a)
    a = 10 ** a

    distances.append(D)

    magnetic_fields.append(B)

    D = D * 9.46 * 10 ** 15  # conversion to meters

    nu = max_freq(B)
    assert nu > 0, f"Maximum emission frequency must be positive, instead got {nu=}."
    nu /= 10**6
    frequencies.append(nu)

    I = complete(B, a, M_s, Mdot, D)
    assert I > 0, f"Radio brightness must be positive, instead got {I=}."

    if IsBurst:
        I = I*(10 ** 1.53)

    intensities.append(I)

    if 30 <= nu <= 180:
        for i in range(6):
            if Freq[i] <= nu <= Freq[i + 1]:
                if I >= (L_EU[i+1] - L_EU[i]) / (Freq[i+1] - Freq[i]) * (nu - Freq[i]) + L_EU[i]:
                    obs = str(exo.name)
                    print(obs, Freq[i], L_EU[i])

    if 120 <= nu < 850 or 1050 < nu <= 1450:
        for i in range(4):
            if uGMRT["Frequencies"][i][0] <= nu <= uGMRT["Frequencies"][i][1] and I > uGMRT["RMS Noise"][i][0]:
                obs = str(exo.name)
                print(obs)

    if exo.name == "tau Boo b":
        obs = exo.name
    if I > 10**(-2):
        obs = exo.name

    labels.append(obs)

arr = np.array(labels)

frequencies = np.array(frequencies)

distances = np.array(distances)
distances = np.reciprocal(distances)
distances *= 10 ** 2.7

intensities = np.array(intensities)

semis = np.array(semis)


df1 = pd.DataFrame({"x": frequencies,
                    "y": intensities,
                    "d": distances,
                    "s": semis,
                    "l": labels})

print(uGMRT)
print(lofar)
# print(lofar.dtypes)

fig0, ax0 = plt.subplots()

im = ax0.scatter(df1.x, df1.y, c=df1.s, s=df1.d, cmap="jet_r")
for i, txt in enumerate(labels):
    if txt:
        ax0.annotate(txt, (df1.x[i], df1.y[i]), fontsize=8)
# lofar.plot(ax=ax0, x="Freq.", y="NL Core", style ="--", linewidth=0.2)
# lofar.plot(ax=ax0, x="Freq.", y="Full EU", style="g--", linewidth=0.5)
ax0.plot(Freq, L_EU, "g--", linewidth=0.5)

for i in range(4):
    x = uGMRT["Frequencies"][i]
    y = uGMRT["RMS Noise"][i]
    plt.plot(x, y, "b--", linewidth=0.5)
    if i == 0:
        ax0.fill_between(x, y, 10 ** 6, color="blue", alpha=0.1, label="uGMRT")
    else:
        ax0.fill_between(x, y, 10 ** 6, color="blue", alpha=0.1)

ax0.fill_between(Freq, L_EU, 10**6, color="green", alpha=0.1, label="LOFAR")
plt.legend()
# ax0.fill_between(Freq, L_NL, 10**6, color="green", alpha=0.1)

fig0.colorbar(im, ax=ax0, label="Distance to Host Star ($\log_{10}{\mathrm{(AU)}}$)")

ax0.axvline(x=10, color="black", linestyle="dashed")

ax0.set_xscale("log")
ax0.set_yscale("log")
# ax0.set_xlim(6, 30)
ax0.set_xlim(left=0.5)
ax0.set_ylim(top=10**1)

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


ax1.hist(intensities, bins=TheBins1, edgecolor="black")
if IsBurst:
    ax1.set_xlabel("Intensity of Burst Emission (Jy)")
    ax1.set_title("Histogram of Burst Emission Intensities")
else:
    ax1.set_xlabel("Intensity of Quiescent Emission (Jy)")
    ax1.set_title("Histogram of Quiescent Emission Intensities")

ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_ylabel("Number of Exoplanets")

fig2, ax2 = plt.subplots()

TheBins2 = np.logspace(-4, 4, 17)

ax2.hist(magnetic_fields, bins=TheBins2, edgecolor="black")
ax2.set_xlabel("Magnetic Field Strength at the Surface (Gauss)")
ax2.set_title("Histogram of the Magnetic Field Strengths")

ax2.set_xscale("log")
ax2.set_ylabel("Number of Exoplanets")

plt.show()

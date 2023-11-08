import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import pandas as pd
from tabulate import tabulate
from adjustText import adjust_text

from radio_module import *

rng = np.random.default_rng()

mpl.use('Qt5Agg')

mpl.rcParams["figure.autolayout"] = True
plt.rcParams['figure.figsize'] = [8, 4]

rc = {"font.family": "times new roman",
      "font.size": 11,
      "mathtext.fontset": "stix"}
plt.rcParams.update(rc)

# LOFAR:
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

# ----------------------------

# uGMRT:
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
# ----------------------------
# MWA

exoplanets = []
IsBurst = 1

names = []
frequencies = []
intensities = []
distances = []
semis = []
magnetic_fields = []
labels = []

y_minerr = []
y_maxerr = []
x_maxerr = []
x_minerr = []

#indices
pl_orbper = 2
pl_orbsmax = 6
pl_bmassj = 14
st_mass = 28
st_age = 32

for i in df.iterrows():
    j = i[1]
    name = j[0]

    T = j[pl_orbper]
    T_b = T + j[pl_orbper+1]
    T_a = T + j[pl_orbper+2]
    if np.isnan(T_b) or np.isnan(T_a) or T_a < 0:
        T_a, T_b = T, T
    T_unc = [T_a, T_b]

    a = j[pl_orbsmax]
    a_b = a + j[pl_orbsmax + 1]
    a_a = a + j[pl_orbsmax + 2]
    if np.isnan(a_b) or np.isnan(a_a) or a_a < 0:
        a_a, a_b = a, a
    a_unc = [a_a, a_b]

    R = j[10]

    if j[18] == "Mass" or j[18] == "Msin(i)/sin(i)":
        M = j[pl_bmassj]
    else:
        M = j[pl_bmassj] * 1.15  # Expected Value of the mass based on projected mass
    M_b = M + j[pl_bmassj+1]
    M_a = M + j[pl_bmassj+2]
    if np.isnan(M_b) or np.isnan(M_a) or M_a < 0:
        M_a, M_b = M, M
    M_unc = [M_a, M_b]

    p = j[19]

    M_s = j[st_mass]
    M_s_b = M_s + j[st_mass+1]
    M_s_a = M_s + j[st_mass+2]
    if np.isnan(M_s_b) or np.isnan(M_s_a) or M_s_a < 0:
        M_s_a, M_s_b = M_s, M_s
    M_s_unc = [M_s_a, M_s_b]

    t = j[st_age]
    t_b = t + j[st_age + 1]
    t_a = t + j[st_age + 2]
    if np.isnan(t_b) or np.isnan(t_a) or t_a < 0:
        t_a, t_b = t, t
    t_unc = [t_a, t_b]

    d = j[36]
    d *= 3.261561

    freqs = []
    intenss = []

    for k in range(1000):
        T = rng.random() * (T_unc[1] - T_unc[0]) + T_unc[0]
        a = rng.random() * (a_unc[1] - a_unc[0]) + a_unc[0]
        M = rng.random() * (M_unc[1] - M_unc[0]) + M_unc[0]
        M_s_i = rng.random() * (M_s_unc[1] - M_s_unc[0]) + M_s_unc[0]
        t = rng.random() * (t_unc[1] - t_unc[0]) + t_unc[0]

        highS_Mdot = t ** (-1.23) * 10 ** 3
        lowS_Mdot = t ** (-0.9) * 10 ** 3
        Mdot = 8.1 * t ** (-1.37)

        B = 1  # Tentative
        sigma = 1  # Jupiter conductivity
        # d *= 3.26156

        exo = Exoplanet(name, a, R, M, p, B, M_s_i, Mdot, d)

        p_c = density(p)
        w_p = rotation(T, a)

        r_c = convective_radius(M, p_c, R)
        mu = magnetic_moment(p_c, w_p, r_c, sigma)
        B = magnetic_field(mu, exo.radius)
        exo.magnetic_field = B
        if B == 0:
            continue

        D = d * 9.46 * 10 ** 15  # conversion to meters

        nu = max_freq(B)
        assert nu > 0, f"Maximum emission frequency must be positive, instead got {nu=}."
        nu /= 10 ** 6
        freqs.append(nu)

        I = complete(B, a, M_s_i, Mdot, D)
        assert I > 0, f"Radio brightness must be positive, instead got {I=}."

        if IsBurst:
            I = I * (10 ** 1.53)

        intenss.append(I)

    y_maxerr.append(max(intenss) - I)
    y_minerr.append(I - min(intenss))
    x_maxerr.append(max(freqs) - nu)
    x_minerr.append(nu - min(freqs))

    nu = np.array(freqs).mean()
    I = np.array(intenss).mean()

    obs = ""

    EXO = Exoplanet(name, a, R, M, p, B, M_s_i, Mdot, d, freq=nu, intensity=I)
    if EXO.magnetic_field != 0:
        exoplanets.append(EXO)

    names.append(name)
    a = np.log10(a)
    semis.append(a)
    distances.append(d)
    magnetic_fields.append(B)
    frequencies.append(nu)
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

    labels.append(obs)

    if exo.name == "tau Boo b":
        plt.hist(intenss)
        # plt.xscale("log")
        # plt.yscale("log")
        plt.show()

arr = np.array(labels)
frequencies = np.array(frequencies)
distances = np.array(distances)
distances = np.reciprocal(distances)
distances *= 10 ** 2.7

intensities = np.array(intensities)

semis = np.array(semis)

y_err = [y_minerr, y_maxerr]
x_err = [x_minerr, x_maxerr]

for i in range(len(labels)):
    if labels[i] == "":
        y_err[0][i], y_err[1][i] = 0, 0
        x_err[0][i], x_err[1][i] = 0, 0

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
ax0.errorbar(df1.x, df1.y,
             yerr=y_err,
             xerr=x_err,
             fmt="None",
             ecolor="black",
             elinewidth=0.5)

# texts = [plt.text(df1.x[i], df1.y[i], labels[i], ha='center', va='center', fontsize=8) for i in range(len(labels)) if labels[i] != ""]
# adjust_text(texts)
for i, txt in enumerate(labels):
    if txt:
        ax0.annotate(txt, xy=(df1.x[i], df1.y[i]), xytext=(1,1), textcoords="offset pixels", fontsize=8)

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

l1 = [frequencies, intensities]
l1 = [arr.tolist() for arr in l1]
l2 = [names]
l2.extend(l1)

table = list(zip(*l2))
file_names = ["names.txt", "freq.txt", "intens.txt"]

table = sorted(table, key=lambda x: x[0].lower())
file_name = file_names[0]
with open(file_name, 'w') as f:
    f.write(tabulate(table))
    f.close()

table = sorted(table, key=lambda x: x[1])
file_name = file_names[1]
with open(file_name, 'w') as f:
    f.write(tabulate(table))
    f.close()

table = sorted(table, key=lambda x: x[2], reverse=True)
file_name = file_names[2]
with open(file_name, 'w') as f:
    f.write(tabulate(table))
    f.close()


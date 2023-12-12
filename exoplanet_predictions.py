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
plt.rcParams['figure.figsize'] = [10, 5]

rc = {"font.family": "times new roman",
      "font.size": 11,
      "mathtext.fontset": "stix"}
plt.rcParams.update(rc)


def retro_noir(ax):
    ax.grid(alpha=0.2)
    ax.tick_params(direction="in", length=7, right=True, top=True, width=1.5)
    ax.spines['top'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['right'].set_linewidth(1.5)


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


# indices
name = 0 #
pl_orbper = 2 #
pl_orbsmax = 6 #
radius = 10 #
pl_bmassj = 14 #
pl_massprov = 18 #
dens = 19 #
st_mass = 28 #
st_age = 32 #
distance = 36

(names, orbs, orb1, orb2, rads, rad1, rad2, smas, smas1, smas2, ms, ms1, ms2,
 massprov, rhos, rho1, rho2, Ms, Ms1, Ms2, ts, t1, t2, ds, d1, d2) = np.genfromtxt("NASA0808.csv",
                                                   usecols=(name, pl_orbper, pl_orbper+1, pl_orbper+2, radius, radius+1, radius+2,
                                                            pl_orbsmax, pl_orbsmax+1, pl_orbsmax+2, pl_bmassj, pl_bmassj+1, pl_bmassj+2,
                                                            pl_massprov, dens, dens+1, dens+2, st_mass, st_mass+1, st_mass+2,
                                                            st_age, st_age+1, st_age+2, distance, distance+1, distance+2),
                                                   skip_header=56, filling_values=0, delimiter=",", unpack=True)
hists = [orbs, smas, ms, Ms, ts, ds]

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

fig, axs = plt.subplots(2, 3, figsize=(12, 5))
plt.rcParams['font.size'] = 8

ax1 = axs[0, 0]
ax2 = axs[0, 1]
ax3 = axs[0, 2]
ax4 = axs[1, 0]
ax5 = axs[1, 1]
ax6 = axs[1, 2]

axes = [ax1, ax2, ax3, ax4, ax5, ax6]

bins = [np.logspace(-1, 6, 15), np.logspace(-3, 3, 13), np.logspace(-3, 2, 11),
        np.logspace(-2, 1, 11), np.logspace(-2, 2, 13), np.logspace(1, 2, 11)]

ax1.hist(orbs, bins=bins[0], edgecolor="black")
ax2.hist(smas, bins=bins[1], edgecolor="black")
ax3.hist(ms,  bins=bins[2], edgecolor="black")
ax4.hist(Ms,  bins=bins[3], edgecolor="black")
ax5.hist(ts,  bins=bins[4], edgecolor="black")
ax6.hist(ds, bins=bins[5], edgecolor="black")

xlabels = ["Orbital Period (Days)", "Semi-major Axis (AU)", f"Planet Mass ($M_j$)", "Star Mass ($M_\odot$)", "Star Age (Gyr)", "Distance (pc)"]
for i in range(len(axes)):
    axes[i].set_xlabel(xlabels[i])
    axes[i].set_xscale("log")
    axes[i].set_yscale("log")

fig.text(0.02, 0.30, 'Bin Count', va='center', rotation='vertical', fontsize=10)
fig.text(0.02, 0.75, 'Bin Count', va='center', rotation='vertical', fontsize=10)

fig.supylabel(' \n', va='center', rotation='vertical', fontsize=11)
fig.suptitle('Distribution of Initial Parameters for the Exoplanet Sample', fontsize=13)

plt.show()

plt.rcParams['font.size'] = 11


for i in df.iterrows():  # The loop that reads exoplanets from NASA file
    j = i[1]
    name = j[0]

    T_i = j[pl_orbper]
    T_s = (j[pl_orbper+1] - j[pl_orbper+2]) / 2
    if np.isnan(T_s):
        T_s = 0

    a_i = j[pl_orbsmax]
    a_s = (j[pl_orbsmax+1] - j[pl_orbsmax+2]) / 2
    if np.isnan(a_s):
        a_s = 0

    R = j[10]

    if j[18] == "Mass" or j[18] == "Msin(i)/sin(i)":
        M_i = j[pl_bmassj]
    else:
        M_i = j[pl_bmassj] * 1.15  # Expected Value of the mass based on projected mass
    M_ss = (j[pl_bmassj + 1] - j[pl_bmassj + 2]) / 2
    if np.isnan(M_ss):
        M_ss = 0

    p = j[19]

    M_s_i = j[st_mass]
    M_s_s = (j[st_mass+1] - j[st_mass+2]) / 2
    if np.isnan(M_s_s):
        M_s_s = 0

    t_i = j[st_age]
    t_s = (j[st_age+1] - j[st_age+2]) / 2
    if np.isnan(t_s):
        t_s = 0

    d = j[36]
    d *= 3.261561

    freqs = []
    intenss = []

    for k in range(1000):  # The loop for Monte Carlo iterations
        T = rng.normal(T_i, T_s)
        while T < 0:
            T = rng.normal(T_i, T_s)

        a = rng.normal(a_i, a_s)
        while a < 0:
            a = rng.normal(a_i, a_s)

        M = rng.normal(M_i, M_ss)
        while M < 0:
            M = rng.normal(M_i, M_ss)

        M_s = rng.normal(M_s_i, M_s_s)
        while M_s < 0:
            M_s = rng.normal(M_i, M_s_s)

        t = rng.normal(t_i, t_s)
        while t < 0:
            t = rng.normal(t_i, t_s)

        highS_Mdot = t ** (-1.23) * 10 ** 3
        lowS_Mdot = t ** (-0.9) * 10 ** 3
        Mdot = 8.1 * t ** (-1.37)

        B = 1  # Tentative
        sigma = 1  # Jupiter conductivity
        # d *= 3.26156

        exo = Exoplanet(name, a, R, M, p, B, M_s, Mdot, d)

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

    nu = np.percentile(freqs, 50)
    I = np.percentile(intenss, 50)

    y_maxerr.append(np.percentile(intenss, 84) - I)
    y_minerr.append(I - np.percentile(intenss, 16))
    x_maxerr.append(np.percentile(freqs, 84) - nu)
    x_minerr.append(nu - np.percentile(freqs, 16))

    obs = ""

    EXO = Exoplanet(name, a, R, M, p, B, M_s, Mdot, d, freq=nu, intensity=I)
    if EXO.magnetic_field != 0:
        exoplanets.append(EXO)

    names.append(EXO.name)
    a = np.log10(EXO.semi_major_axis)
    semis.append(a)
    distances.append(EXO.distance)
    magnetic_fields.append(EXO.magnetic_field)
    frequencies.append(EXO.freq)
    intensities.append(EXO.intensity)

    if 30 <= nu <= 180:
        for m in range(6):
            if Freq[m] <= nu <= Freq[m + 1]:
                if I >= (L_EU[m + 1] - L_EU[m]) / (Freq[m + 1] - Freq[m]) * (nu - Freq[m]) + L_EU[m]:
                    obs = str(EXO.name)
                    print(obs)

    if 120 <= nu < 850 or 1050 < nu <= 1450:
        for m in range(4):
            if uGMRT["Frequencies"][m][0] <= nu <= uGMRT["Frequencies"][m][1] and I > uGMRT["RMS Noise"][m][0]:
                obs = str(EXO.name)
                print(obs)

    if EXO.name == "tau Boo b":  # Special interest
        obs = EXO.name

    labels.append(obs)

    selected = "tau Boo b"
    if EXO.name == selected:
        plt.subplot(1, 2, 1)
        plt.hist(intenss, edgecolor="black")
        # plt.xscale("log")
        # plt.yscale("log")
        plt.title(f"Distribution of Emission Intensity for {selected}")
        plt.xlabel("Intensity (Jy)")
        plt.ylabel("Bin Count")

        plt.subplot(1, 2, 2)
        plt.hist(freqs, edgecolor="black")
        # plt.xscale("log")
        # plt.yscale("log")
        plt.title(f"Distribution of Emission Frequency for {selected}")
        plt.xlabel("Frequency (MHz)")
        plt.ylabel("Bin Count")
        # plt.show()
        # plt.close()


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

# print(uGMRT)
# print(lofar)
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
        ax0.annotate(txt, xy=(df1.x[i], df1.y[i]), xytext=(1, 1), textcoords="offset pixels", fontsize=8)

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

fig1, axs = plt.subplots(1, 2, sharey="row", figsize=[10, 5])

ax1, ax2 = axs[0], axs[1]

ax1.hist(intensities, bins=TheBins1, edgecolor="black")
if IsBurst:
    ax1.set_xlabel("Intensity of Burst Emission (Jy)")
    ax1.set_title("Histogram of Burst Emission Intensities")
else:
    ax1.set_xlabel("Intensity of Quiescent Emission (Jy)")
    ax1.set_title("Histogram of Quiescent Emission Intensities")

ax1.set_xscale("log")
ax1.set_yscale("log")

TheBins2 = np.logspace(-4, 4, 17)

ax2.hist(magnetic_fields, bins=TheBins2, edgecolor="black")
ax2.set_xlabel("Magnetic Field Strength at the Surface (Gauss)")
ax2.set_title("Histogram of the Magnetic Field Strengths")

ax2.set_xscale("log")
fig1.supylabel("Number of Exoplanets")

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

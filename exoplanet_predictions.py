import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
from tabulate import tabulate
from adjustText import adjust_text
from radio_module import *
from rotation_script import *
import smplotlib
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D
from tqdm import tqdm
from astropy.coordinates import SkyCoord
import astropy.units as u

plt.rcParams['figure.figsize'] = [10, 5]

rng = np.random.default_rng()

# LOFAR:
lofar = pd.read_csv(
    "sensitivities.csv")  # Obtained from van Haarlem et al. (2017), 8h integration time, 4.66MHz effective bandwidth
lofar.columns = lofar.iloc[0]
lofar = lofar[1:]
lofar = lofar.drop([1, 2, 3, 11, 12, 13])
lofar = lofar.reset_index(drop=True)
lofar = lofar.apply(pd.to_numeric, errors="ignore")
lofar["NL Core"] = lofar["NL Core"].multiply(10 ** (-3)) * 5  # 5 sigma sensitivity in Jy
lofar["Full EU"] = lofar["Full EU"].multiply(10 ** (-3)) * 5  # 5 sigma sensitivity in Jy

L_NL = lofar["NL Core"]
L_EU = lofar["Full EU"]
Freq = lofar["Freq."]

L_EU_1 = L_EU[:4]
L_EU_2 = L_EU[4:]
Freq_1 = Freq[:4]
Freq_2 = Freq[4:]

# ----------------------------
# NenuFAR
NenuNoise = np.array([130.0, 9.0])  # MJy (10 MHz, 1h, website)
NenuFreq = np.array([15.0, 85.0])  # MHz (website)
NenuNoise *= np.sqrt(1 / 8) * 10 ** (-3) * 5  # 5 sigma sensitivity in Jy
nenu_data = np.load("nenufar.npz")
NenuFreq = nenu_data["nenu_freqs"]
NenuNoise = nenu_data["nenu_noise"] * 5

# ----------------------------
# uGMRT:
d = {"Bands": ["Band 1", "Band 2", "Band 3", "Band 4"],
     "Frequencies": [[120, 250], [250, 500], [550, 850], [1050, 1450]],  # MHz
     "RMS Noise": [np.array([190, 190]), np.array([50, 50]), np.array([40, 40]), np.array([45, 45])]
     # microJy, 10min integration time, 100MHz Bandwidth
     }

uGMRT = pd.DataFrame(data=d)
integration_time = 8 * 60  # minutes
bandwidth = 100  # MHZ
uGMRT["RMS Noise"] = uGMRT["RMS Noise"] * (np.sqrt((100 * 10) / (bandwidth * integration_time))) * 10 ** (
    -6) * 5  # 5 sigma sensitivity in Jy

# ----------------------------
# MWA

d_mwa = {"Frequencies": [[72.30, 103.04], [103.04, 133.76], [138.88, 169.60], [169.60, 200.32], [200.32, 231.04]],
         "RMS Noise": [np.array([24, 24]), np.array([14, 14]), np.array([7, 7]), np.array([5, 5]),
                       np.array([5, 5])]}  # mJy, 2min integration time, 40kHz bandwidth

MWA = pd.DataFrame(data=d_mwa)
integration_time = 8 * 60  # minutes
MWA["RMS Noise"] *= np.sqrt((2 / integration_time)) * 10 ** (-3) * 5

# Retrieve Data
filename = "NASA1904.csv"
df = pd.read_csv(filename, comment="#")
windfile = "new_wind_info" + filename[4:-4] + ".txt"

headers_to_find = ["pl_name", "pl_orbper", "pl_orbsmax", "pl_radj", "pl_bmassj", "pl_bmassprov", "pl_dens", "st_rad",
                   "st_mass", "st_age", "sy_dist"]
indices = [df.columns.get_loc(header) for header in headers_to_find]
pl_name, pl_orbper, pl_orbsmax, radius, pl_bmassj, pl_bmassprov, dens, st_rad, st_mass, st_age, distance = indices

wind_temperatures, wind_speeds = np.genfromtxt(windfile, usecols=(1, 2), skip_header=1, delimiter="\t", unpack=True)

# Plot initial distributions
orbs = df["pl_orbper"].to_numpy()
smas = df["pl_orbsmax"].to_numpy()
ms = df["pl_bmassj"].to_numpy()
Ms = df["st_mass"].to_numpy()
ts = df["st_age"].to_numpy()
ds = df["sy_dist"].to_numpy()

hists = [orbs, smas, ms, Ms, ts, ds]

exoplanets = []
IsBurst = 1

names = []
frequencies = []
intensities_mag = []
intensities_kin = []
intensities_both = []
distances = []
semis = []
magnetic_fields = []
labels_mag = []
labels_kin = []
labels_both = []
outliers = []

y_mag_minerr = []
y_mag_maxerr = []
y_kin_minerr = []
y_kin_maxerr = []
y_both_maxerr = []
y_both_minerr = []
x_maxerr = []
x_minerr = []

fig, axs = plt.subplots(2, 3, figsize=(10, 5))
plt.rcParams['font.size'] = 12

ax1 = axs[0, 0]
ax2 = axs[0, 1]
ax3 = axs[0, 2]
ax4 = axs[1, 0]
ax5 = axs[1, 1]
ax6 = axs[1, 2]

axes = [ax1, ax2, ax3, ax4, ax5, ax6]

bins = [np.logspace(-1, 5, 13), np.logspace(-2.5, 1.5, 13), np.logspace(-3, 1.5, 9),
        np.logspace(-1, 1, 9), np.logspace(-2, 1.5, 15), np.logspace(0.5, 3.5, 10)]

ax1.hist(orbs, bins=bins[0], edgecolor="black", color="xkcd:ocean")
ax2.hist(smas, bins=bins[1], edgecolor="black", color="xkcd:ocean")
ax3.hist(ms, bins=bins[2], edgecolor="black", color="xkcd:ocean")
ax4.hist(Ms, bins=bins[3], edgecolor="black", color="xkcd:ocean")
ax5.hist(ts, bins=bins[4], edgecolor="black", color="xkcd:ocean")
ax6.hist(ds, bins=bins[5], edgecolor="black", color="xkcd:ocean")

xlabels = ["Orbital Period (Days)", "Semi-major Axis (AU)", f"Planet Mass ($M_j$)", "Star Mass ($M_\odot$)",
           "Star Age (Gyr)", "Distance (pc)"]
for i in range(len(axes)):
    axes[i].set_xlabel(xlabels[i])
    axes[i].set_xscale("log")
    axes[i].set_yscale("log")

fig.text(0, 0.30, 'Bin Count', va='center', rotation='vertical', fontsize=12)
fig.text(0, 0.75, 'Bin Count', va='center', rotation='vertical', fontsize=12)

fig.tight_layout()

fig.supylabel(' \n', va='center', rotation='vertical', fontsize=11)
fig.savefig("dist.pdf")
plt.show()
plt.close(fig)

plt.rcParams['font.size'] = 10

detectables_mag = []
detectables_kin = []
detectables_both = []
lofar_obs = []
nenufar_obs = []
mwa_obs = []
ugmrt_obs = []
x_low, x_up, y_low, y_up = [], [], [], []
outliers = set()
insiders = set()

# Calculate Frequencies and intensities
for i, j in tqdm(df.iterrows(), total=len(df)):  # The loop that reads exoplanets from NASA file

    name = j[pl_name]

    T_i = j[pl_orbper]
    T_s = (j[pl_orbper + 1] - j[pl_orbper + 2]) / 2
    if np.isnan(T_s):
        T_s = T_i / 5

    a_i = j[pl_orbsmax]
    a_s = (j[pl_orbsmax + 1] - j[pl_orbsmax + 2]) / 2
    if np.isnan(a_s):
        a_s = a_i / 5

    R_i = j[radius]
    R_s = (j[radius + 1] - j[radius + 2]) / 2
    if np.isnan(R_s):
        R_s = R_i / 5

    if j[pl_bmassprov] == "Mass" or j[pl_bmassprov] == "Msin(i)/sin(i)":
        M_i = j[pl_bmassj]
    else:
        M_i = j[pl_bmassj] * 1.15  # Expected Value of the mass based on projected mass
    M_ss = (j[pl_bmassj + 1] - j[pl_bmassj + 2]) / 2
    if np.isnan(M_ss):
        M_ss = M_i / 5

    p_i = j[dens]
    p_s = (j[dens + 1] - j[dens + 2]) / 2
    if np.isnan(p_s):
        p_s = p_i / 5

    M_s_i = j[st_mass]
    M_s_s = (j[st_mass + 1] - j[st_mass + 2]) / 2
    if np.isnan(M_s_s):
        M_s_s = M_s_i / 5

    t_i = j[st_age]
    t_s = (j[st_age + 1] - j[st_age + 2]) / 2
    if np.isnan(t_s):
        t_s = t_i / 5

    Rs_i = j[st_rad]
    Rs_s = (j[st_rad + 1] - j[st_rad + 2]) / 2
    if np.isnan(Rs_s):
        Rs_s = Rs_i / 5

    T_wind = wind_temperatures[i]
    v_wind = wind_speeds[i]

    d = j[distance]
    d *= 3.261561

    freqs = []
    intens_mag = []
    intens_kin = []
    intens_both = []
    n_ps = []

    high_var = ["AU Mic c", "V1298 Tau d", "V1298 Tau b", "V1298 Tau e", "V1298 Tau c"]

    for k in range(10000):  # The loop for Monte Carlo iterations
        T = rng.normal(T_i, T_s)
        while T < 0:
            T = rng.normal(T_i, T_s)

        a = rng.normal(a_i, a_s)
        while a < 0:
            a = rng.normal(a_i, a_s)

        R = rng.normal(R_i, R_s)
        while R < 0:
            R = rng.normal(R_i, R_s)

        M = rng.normal(M_i, M_ss)
        while M < 0:
            M = rng.normal(M_i, M_ss)

        p = rng.normal(p_i, p_s)
        while p < 0:
            p = rng.normal(p_i, p_s)

        M_s = rng.normal(M_s_i, M_s_s)
        while M_s < 0:
            M_s = rng.normal(M_i, M_s_s)

        t = rng.normal(t_i, t_s)
        while t < 0:
            t = rng.normal(t_i, t_s)

        Rs = rng.normal(Rs_i, Rs_s)
        while Rs < 0:
            Rs = rng.normal(Rs_i, Rs_s)

        b0_exponent = rng.normal(-0.655, 0.045)  # Vidotto 2014

        flux_exponent = rng.normal(-1.74, 0.34)  # Ayres 1997 in Lynch 2018
        loss_exponent = rng.normal(0.79, 0.17)  # Alvarado-Gomez 2016 in Lynch 2018

        L = moment_sampler()

        highS_Mdot = t ** (-1.23) * 10 ** 3
        lowS_Mdot = t ** (-0.9) * 10 ** 3
        Mdot = mass_loss(t, flux_exponent, loss_exponent)

        sigma = 1  # Jupiter conductivity

        p_c = density(p)
        w_p = rotation(T, a, L, M, R)

        r_c = convective_radius(M, p_c, R)
        mu = magnetic_moment(p_c, w_p, r_c, sigma)
        B = magnetic_field(mu, R)

        if B == 0:
            continue

        D = d * 9.46 * 10 ** 15  # conversion to meters

        B_perp = imf_perp_complete(M_s, a, Rs, t, v_wind, b0_exponent)
        v_k = keplerian(M_s, a)
        veff = v_eff(v_wind, v_k)
        B_star = imf_complete(M_s, a, v_wind, t, Rs, b0_exponent)
        n = number_density(Mdot, veff, a)
        R_m = Rm(B, R, n, T_wind, v_wind, B_star)

        n_p = 8.98 * np.sqrt(n) * 10 ** (-3)

        nu = max_freq(B)
        assert nu > 0, f"Maximum emission frequency must be positive, instead got {nu=}."

        nu /= 10 ** 6
        n_p /= 10 ** 6

        flag = True

        if freq_condition(nu, n) or name in high_var:
            flag = False
            freqs.append(nu)
            n_ps.append(n_p)
            nu *= 10 ** 6

            I = complete(B, a, M_s, Mdot, D)
            assert I > 0, f"Radio brightness must be positive, instead got {I=}."

            P_in_mag = P_input_mag(B_perp, veff * 10 ** 3, R_m * 7 * 10 ** 8, n)
            P_in_kin = P_input_kin(B_perp, veff * 10 ** 3, R_m * 7 * 10 ** 8, n)
            P_rad_mag = radio_power(P_in_mag, 0, nu, D)
            P_rad_kin = radio_power(0, P_in_kin, nu, D)
            P_rad_both = radio_power(P_in_mag, P_in_kin, nu, D, both=True)

            if IsBurst:
                factor = 10  # 10**1.53 (Ashtari)
                I = I * factor
                P_rad_mag *= factor
                P_rad_kin *= factor
                P_rad_both *= factor

            intens_mag.append(P_rad_mag)
            intens_kin.append(P_rad_kin)
            intens_both.append(P_rad_both)

    if flag:
        continue

    freqs = np.array(freqs)
    n_p = np.percentile(n_ps, 50)
    nu = np.percentile(freqs, 50)
    I_mag = np.percentile(intens_mag, 50)
    I_kin = np.percentile(intens_kin, 50)
    I_both = np.percentile(intens_both, 50)

    y_mag_maxerr.append(np.percentile(intens_mag, 84) - I_mag)
    y_mag_minerr.append(I_mag - np.percentile(intens_mag, 16))
    y_kin_maxerr.append(np.percentile(intens_kin, 84) - I_kin)
    y_kin_minerr.append(I_kin - np.percentile(intens_kin, 16))
    y_both_maxerr.append(np.percentile(intens_both, 84) - I_both)
    y_both_minerr.append(I_both - np.percentile(intens_both, 16))

    x_maxerr.append(np.percentile(freqs, 84) - nu)
    x_minerr.append(nu - np.percentile(freqs, 16))

    x_low.append((np.percentile(freqs, 16)) / nu)
    x_up.append((np.percentile(freqs, 84)) / nu)

    y_low.append((np.percentile(intens_both, 16)) / I_both)
    y_up.append((np.percentile(intens_both, 84)) / I_both)

    x_err_avg = (-np.log10((np.percentile(freqs, 16)) / nu) + np.log10((np.percentile(freqs, 84)) / nu)) / 2
    y_err_avg = (-np.log10((np.percentile(intens_both, 16)) / I_both) + np.log10(
        (np.percentile(intens_both, 84)) / I_both)) / 2

    obs_mag = ""
    obs_kin = ""
    obs_both = ""
    out = ""

    EXO = Exoplanet(name, a, R, M, p, B, M_s, Mdot, d, freq=nu, intensity_mag=I_mag, intensity_kin=I_kin,
                    intensity_both=I_both)
    if EXO.magnetic_field != 0:
        exoplanets.append(EXO)

    names.append(EXO.name)
    a = np.log10(EXO.semi_major_axis)
    semis.append(a)
    distances.append(EXO.distance)
    magnetic_fields.append(EXO.magnetic_field)
    frequencies.append(EXO.freq)
    intensities_mag.append(EXO.intensity_mag)
    intensities_kin.append(EXO.intensity_kin)
    intensities_both.append(EXO.intensity_both)

    observable_mag = False
    observable_kin = False
    observable_both = False
    observable_flag = False

    if x_err_avg < 0.7 and nu > 10:
        frac_lim = 0.3

        for m in range(3):
            x = freqs[(Freq_1[m] <= freqs) & (freqs <= Freq_1[m + 1])]
            frac = len(x) / len(freqs)

            if frac > frac_lim:
                # if I_mag >= (L_EU_1[m + 1] - L_EU_1[m]) / (Freq_1[m + 1] - Freq_1[m]) * (nu - Freq_1[m]) + L_EU_1[m]:
                if I_mag >= (L_EU_1[m] + L_EU_1[m + 1]) / 2:
                    obs_mag = str(EXO.name)
                    # print(obs_mag, m, I_mag)
                    observable_mag = True
                # if I_kin >= (L_EU_1[m + 1] - L_EU_1[m]) / (Freq_1[m + 1] - Freq_1[m]) * (nu - Freq_1[m]) + L_EU_1[m]:
                if I_kin >= (L_EU_1[m] + L_EU_1[m + 1]) / 2:
                    obs_kin = str(EXO.name)
                    # print(obs_kin, m, I_kin)
                    observable_kin = True
                # if I_both >= (L_EU_1[m + 1] - L_EU_1[m]) / (Freq_1[m + 1] - Freq_1[m]) * (nu - Freq_1[m]) + L_EU_1[m]:
                if I_both >= (L_EU_1[m] + L_EU_1[m + 1]) / 2:
                    obs_both = str(EXO.name)
                    # print(obs_both, m, I_both)
                    observable_both = True
                    lofar_obs.append(obs_both)
                    if not Freq_1[0] <= nu <= Freq_1[3]:
                        out = str(EXO.name)
                        # print(f"{out} Outlier for LOFAR LBA {m=}")
                        outliers.add(obs_both)
                    else:
                        insiders.add(obs_both)
                observable_flag = observable_both or observable_mag or observable_kin

        for m in range(4, 6):
            x = freqs[(Freq_2[m] <= freqs) & (freqs <= Freq_2[m + 1])]
            frac = len(x) / len(freqs)

            if frac > frac_lim:
                # if I_mag >= (L_EU_2[m + 1] - L_EU_2[m]) / (Freq_2[m + 1] - Freq_2[m]) * (nu - Freq_2[m]) + L_EU_2[m]:
                if I_mag >= (L_EU_2[m] + L_EU_2[m + 1]) / 2:
                    obs_mag = str(EXO.name)
                    # print(obs_mag)
                    observable_mag = True

                # if I_kin >= (L_EU_2[m + 1] - L_EU_2[m]) / (Freq_2[m + 1] - Freq_2[m]) * (nu - Freq_2[m]) + L_EU_2[m]:
                if I_kin >= (L_EU_2[m] + L_EU_2[m + 1]) / 2:
                    obs_kin = str(EXO.name)
                    # print(obs_kin)
                    observable_kin = True
                # if I_both >= (L_EU_2[m + 1] - L_EU_2[m]) / (Freq_2[m + 1] - Freq_2[m]) * (nu - Freq_2[m]) + L_EU_2[m]:

                if I_both >= (L_EU_2[m] + L_EU_2[m + 1]) / 2:
                    obs_both = str(EXO.name)
                    # print(obs_both)
                    observable_both = True
                    lofar_obs.append(obs_both)
                    if not Freq_2[4] <= nu <= Freq_2[6]:
                        out = str(EXO.name)
                        # print(f"{out} Outlier for LOFAR HBA {m=} ")
                        outliers.add(obs_both)
                    else:
                        insiders.add(obs_both)
                observable_flag = observable_both or observable_mag or observable_kin

        # NenuFAR
        x = freqs[(NenuFreq[0] <= freqs) & (freqs <= NenuFreq[-1])]
        frac = len(x) / len(freqs)
        for i in range(len(NenuFreq) - 1):

            if NenuFreq[i] <= nu <= NenuFreq[i + 1]:
                if I_mag >= NenuNoise[i] + (NenuNoise[i + 1] - NenuNoise[i]) / (NenuFreq[i + 1] - NenuFreq[i]) * (
                        nu - NenuFreq[i]):
                    obs_mag = str(EXO.name)
                    # print(obs_mag)
                    observable_mag = True

                if I_kin >= NenuNoise[i] + (NenuNoise[i + 1] - NenuNoise[i]) / (NenuFreq[i + 1] - NenuFreq[i]) * (
                        nu - NenuFreq[i]):
                    obs_kin = str(EXO.name)
                    # print(obs_kin)
                    observable_kin = True

                if I_both >= NenuNoise[i] + (NenuNoise[i + 1] - NenuNoise[i]) / (NenuFreq[i + 1] - NenuFreq[i]) * (
                        nu - NenuFreq[i]):
                    obs_both = str(EXO.name)
                    # print(obs_both)
                    observable_both = True
                    nenufar_obs.append(obs_both)
                    insiders.add(obs_both)

            elif frac > frac_lim and I_both > NenuNoise[0] * 2 / 3:
                obs_both = str(EXO.name)
                # print(obs_both)
                observable_both = True
                nenufar_obs.append(obs_both)
                out = str(EXO.name)
                # print(f"{out} Outlier for NenuFAR")
                outliers.add(obs_both)
            # observable_flag = observable_both or observable_mag or observable_kin

        for m in range(5):
            x = freqs[(MWA["Frequencies"][m][0] <= freqs) & (freqs <= MWA["Frequencies"][m][1])]
            frac = len(x) / len(freqs)

            if frac > frac_lim and I_mag > MWA["RMS Noise"][m][0]:
                obs_mag = str(EXO.name)
                # print(obs_mag)
                observable_mag = True

            if frac > frac_lim and I_kin > MWA["RMS Noise"][m][0]:
                obs_kin = str(EXO.name)
                # print(obs_kin)
                observable_kin = True

            if frac > frac_lim and I_both > MWA["RMS Noise"][m][0]:
                obs_both = str(EXO.name)
                # print(obs_both)
                observable_both = True
                mwa_obs.append(obs_both)
                if not MWA["Frequencies"][0][0] <= nu <= MWA["Frequencies"][4][1]:
                    out = str(EXO.name)
                    # print(f"{out} Outlier for MWA {m=} ")
                    outliers.add(obs_both)
                else:
                    insiders.add(obs_both)
                break

        for m in range(4):

            x = freqs[(uGMRT["Frequencies"][m][0] <= freqs) & (freqs <= uGMRT["Frequencies"][m][1])]
            frac = len(x) / len(freqs)

            if frac > frac_lim and I_mag > uGMRT["RMS Noise"][m][0]:
                obs_mag = str(EXO.name)
                # print(obs_mag)
                observable_mag = True

            if frac > frac_lim and I_kin > uGMRT["RMS Noise"][m][0]:
                obs_kin = str(EXO.name)
                # print(obs_kin)
                observable_kin = True

            if frac > frac_lim and I_both > uGMRT["RMS Noise"][m][0]:
                obs_both = str(EXO.name)
                # print(obs_both)
                observable_both = True
                ugmrt_obs.append(obs_both)
                if not uGMRT["Frequencies"][0][0] <= nu <= uGMRT["Frequencies"][3][1]:
                    out = str(EXO.name)
                    # print(f"{out} Outlier for uGMRT {m=}")
                    outliers.add(obs_both)
                else:
                    insiders.add(obs_both)
                break

    # if str(EXO.name) in high_var and not observable_both:
    #     obs_both = str(EXO.name)
    #     detectables_both.append(EXO)

    if EXO.name == "tau Boo b":  # Special interest
        obs = EXO.name

    if observable_both:
        detectables_both.append(EXO)
    if observable_kin:
        detectables_kin.append(EXO)
    if observable_mag:
        detectables_mag.append(EXO)

    labels_mag.append(obs_mag)
    labels_kin.append(obs_kin)
    labels_both.append(obs_both)
    # outliers.append(out)

    labels = [labels_mag, labels_kin, labels_both]

    selected = ["tau Boo b", "AU Mic c", "V1298 Tau d", "V1298 Tau b", "V1298 Tau e", "V1298 Tau c"]
    if EXO.name in selected:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(4, 6.4))
        ax1.hist(np.log10(intens_both), edgecolor="black")
        ax1.axvline(np.log10(I_both), linestyle="--", color="crimson", label="Median")
        ax1.axvline(np.log10(np.percentile(intens_both, 16)), color="xkcd:deep green", linestyle="--",
                    label="16th & 84th\npercentiles")
        ax1.axvline(np.log10(np.percentile(intens_both, 84)), color="xkcd:deep green", linestyle="--")
        # ax1.set_title(f"{EXO.name}")
        ax1.set_xlabel("log$_{10}$(Flux Density [Jy])")
        ax1.set_ylabel("Bin Count")

        ax2.hist(np.log10(freqs), histtype="step")
        if n_p > min(freqs):
            ax2.axvline(np.log10(n_p), linestyle="dotted", label="Local Plasma\nFrequency")
        ax2.axvline(np.log10(nu), linestyle="--", color="firebrick")
        if min(freqs) < 1:
            ax2.axvline(1, linestyle="--", label="Ionosphere\ncutoff")
        ax2.axvline(np.log10(np.percentile(freqs, 16)), color="xkcd:deep green", linestyle="--")
        ax2.axvline(np.log10(np.percentile(freqs, 84)), color="xkcd:deep green", linestyle="--")
        # ax2.set_title(f"{EXO.name}")
        ax2.set_xlabel("log$_{10}$(Frequency [MHz])")
        ax2.set_ylabel("Bin Count")

        fig.suptitle(f"{EXO.name}")
        lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
        lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
        fig.legend(lines, labels, loc="upper right", bbox_to_anchor=(1, 1), ncol=int(len(lines_labels) / 2),
                   frameon=True, shadow=True)
        plt.savefig(f"Distribution Plots/{EXO.name}.pdf")

real_outliers = set(outliers) - set(insiders)
lofar_obs, nenufar_obs, mwa_obs, ugmrt_obs = np.array(lofar_obs), np.array(nenufar_obs), np.array(mwa_obs), np.array(
    ugmrt_obs)
frequencies = np.array(frequencies)
distances = np.array(distances)
distances = np.reciprocal(distances)
distances *= 10 ** 2.7

intensities_mag = np.array(intensities_mag)
intensities_kin = np.array(intensities_kin)
intensities_both = np.array(intensities_both)
intensities = [intensities_mag, intensities_kin, intensities_both]

semis = np.array(semis)
cond = np.where(semis < 2)

y_mag_minerr = np.array(y_mag_minerr)
y_mag_maxerr = np.array(y_mag_maxerr)
y_kin_minerr = np.array(y_kin_minerr)
y_kin_maxerr = np.array(y_kin_maxerr)
y_both_maxerr = np.array(y_both_maxerr)
y_both_minerr = np.array(y_both_minerr)
x_minerr = np.array(x_minerr)
x_maxerr = np.array(x_maxerr)

y_mag_err = [y_mag_minerr, y_mag_maxerr]
y_kin_err = [y_kin_minerr, y_kin_maxerr]
y_both_err = [y_both_minerr, y_both_maxerr]
y_err = [y_mag_err, y_kin_err, y_both_err]
x_err = [x_minerr, x_maxerr]

x_low_avg = np.exp(np.mean(np.log(np.array(x_low))))
x_up_avg = np.exp(np.mean(np.log(np.array(x_up))))
y_low_avg = np.exp(np.mean(np.log(np.array(y_low))))
y_up_avg = np.exp(np.mean(np.log(np.array(y_up))))
average_errors = [x_low_avg, x_up_avg, y_low_avg, y_up_avg]

np.savez("all.npz", y_mag_minerr=y_mag_minerr, y_mag_maxerr=y_mag_maxerr,
         y_kin_maxerr=y_mag_maxerr, y_kin_minerr=y_mag_minerr, y_both_minerr=y_both_minerr,
         y_both_maxerr=y_both_maxerr, x_minerr=x_minerr, x_maxerr=x_maxerr, average_errors=np.array(average_errors),
         detectables_both=np.array(detectables_both), intensities=np.array(intensities),
         magnetic_fields=np.array(magnetic_fields), real_outliers=np.array(list(real_outliers)))


df1 = pd.DataFrame({"Names": names,
                    "x": frequencies,
                    "y_mag": intensities_mag,
                    "y_kin": intensities_kin,
                    "y_both": intensities_both,
                    "d": distances,
                    "s": semis,
                    "l_mag": labels_mag,
                    "l_kin": labels_kin,
                    "l_both": labels_both,
                    "y_err_min": y_both_minerr,
                    "y_err_max": y_both_maxerr,
                    "x_err_min": x_minerr,
                    "x_err_max": x_maxerr})

df1['Color'] = df1['Names'].apply(lambda x: 'gray' if x in real_outliers else 'black').astype(str)
df1.to_csv("df1.csv", index=False)

det_data = [[exo.name, exo.freq, exo.intensity_both] for exo in detectables_both]
df_det = pd.DataFrame(det_data[0:], columns=["Name", "Freq", "Flux"])

for i in range(len(intensities)):
    l1 = [frequencies, intensities[i]]
    l1 = [arr.tolist() for arr in l1]
    l2 = [names]
    l2.extend(l1)

    if i == 0:
        file_names = ["names_mag.txt", "freq_mag.txt", "intens_mag.txt", "detectables_mag.csv"]
        detectables = detectables_mag
        detectable_data = pd.DataFrame([[exo.name, exo.freq, exo.intensity_mag] for exo in detectables],
                                       columns=["Name", "Freq", "Flux"])

    elif i == 1:
        file_names = ["names_kin.txt", "freq_kin.txt", "intens_kin.txt", "detectables_kin.csv"]
        detectables = detectables_kin
        detectable_data = pd.DataFrame([[exo.name, exo.freq, exo.intensity_kin] for exo in detectables],
                                       columns=["Name", "Freq", "Flux"])

    else:
        file_names = ["names_both.txt", "freq_both.txt", "intens_both.txt", "detectables_both.csv"]
        detectables = detectables_both
        detectable_data = pd.DataFrame([[exo.name, exo.freq, exo.intensity_both] for exo in detectables],
                                       columns=["Name", "Freq", "Flux"])

    # with open(file_names[3], "w") as fn:
    #     fn.write(tabulate(detectable_data))
    #     fn.close()

    detectable_data["Flux"] *= 10 ** 3


    def freq_format(x):
        return '{:.2f}'.format(x)


    def flux_format(x):
        return '{:.3f}'.format(x)


    detectable_data["Freq"] = detectable_data["Freq"].map(freq_format)
    detectable_data["Flux"] = detectable_data["Flux"].map(flux_format)
    detectable_data = detectable_data.sort_values(by="Name")
    detectable_data.columns = ["Name", "Freq(MHz)", "Flux(mJy)"]
    detectable_data.to_csv(file_names[3], index=False)

    table = list(zip(*l2))

    table = sorted(table, key=lambda x: x[0].lower())
    file_name = f"Old Result Tables/{file_names[0]}"
    with open(file_name, 'w') as f:
        f.write(tabulate(table))
        f.close()

    table = sorted(table, key=lambda x: x[1])
    file_name = f"Old Result Tables/{file_names[1]}"
    with open(file_name, 'w') as f:
        f.write(tabulate(table))
        f.close()

    table = sorted(table, key=lambda x: x[2], reverse=True)
    file_name = f"Old Result Tables/{file_names[2]}"
    with open(file_name, 'w') as f:
        f.write(tabulate(table))
        f.close()

final_df = df1[["Names", "x", "y_both", "x_err_min", "x_err_max", "y_err_min", "y_err_max"]]
# Step 1: Merge the two dataframes on the name column
merged = final_df.merge(df, left_on='Names', right_on='pl_name', how='left')

# Step 2: Define the formatting function
def format_coords(row):
    coord = SkyCoord(ra=row['ra']*u.deg, dec=row['dec']*u.deg)
    ra_str = coord.ra.to_string(unit=u.hour, sep=':', precision=0, pad=True)
    dec_str = coord.dec.to_string(unit=u.deg, sep=':', precision=0, alwayssign=True, pad=True)
    return pd.Series({'RA_formatted': ra_str, 'DEC_formatted': dec_str})

# Step 3: Apply formatting to merged dataframe
formatted_coords = merged.apply(format_coords, axis=1)

# Step 4: Add formatted columns to final_df
final_df[['RA', 'DEC']] = formatted_coords
final_df.columns = ["Name", "Max. Frequency [MHz]", "Max. Flux Density [Jy]", "Max. Frequency Lower Uncertainty [MHz]",
                    "Max. Frequency Upper Uncertainty [MHz]", "Max. Flux Density Lower Uncertainty [Jy]",
                    "Max. Flux Density Upper Uncertainty [Jy]", "RA (J2000)", "DEC (J2000)"]
# List of columns, moving "RA (J2000)" and "DEC (J2000)" right after "Name"
cols = final_df.columns.tolist()
cols.remove("RA (J2000)")
cols.remove("DEC (J2000)")

# Insert coordinate columns after "Name"
name_index = cols.index("Name")
cols[name_index+1:name_index+1] = ["RA (J2000)", "DEC (J2000)"]
final_df = final_df[cols]

final_df = final_df.sort_values(by="Max. Flux Density [Jy]", ascending=False)

final_df.to_csv("Output Tables/all.csv", index=False)
end0 = final_df[(final_df["Max. Frequency [MHz]"] > 1e-1) & (final_df["Max. Frequency [MHz]"] < 1)]
end1 = final_df[(final_df["Max. Frequency [MHz]"] > 1) & (final_df["Max. Frequency [MHz]"] < 1e1)]
end2 = final_df[(final_df["Max. Frequency [MHz]"] > 1e1) & (final_df["Max. Frequency [MHz]"] < 1e2)]
end3 = final_df[(final_df["Max. Frequency [MHz]"] > 1e2) & (final_df["Max. Frequency [MHz]"] < 1e3)]
end0.to_csv("Output Tables/0.1-1 MHz.csv", index=False)
end1.to_csv("Output Tables/1-10 MHz.csv", index=False)
end2.to_csv("Output Tables/10-100 MHz.csv", index=False)
end3.to_csv("Output Tables/100-1000 MHz.csv", index=False)

np.savez("observables", lofar=lofar_obs, nenufar=nenufar_obs, mwa=mwa_obs, ugmrt=ugmrt_obs)

with open("result_stats.txt", "w") as fl:
    fl.write(f"Total: {len(df)}\n"
             f"Can reach: {len(frequencies)}\n"
             f"Can reach (%): {len(frequencies) / len(df) * 100:.1f}\n"
             f"Can penetrate ionosphere:  {len(frequencies[np.where(frequencies > 10)])}\n"
             f"Can penetrate ionosphere in can reach (%): {len(frequencies[np.where(frequencies > 10)]) / len(frequencies) * 100:.1f}\n"
             f"Not close-in but reach: {len(frequencies[np.where((frequencies > 10) & (semis > -1))])}\n")
    fl.close()

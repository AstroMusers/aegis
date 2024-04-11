import numpy as np
import pandas as pd
import math
from tabulate import tabulate
from adjustText import adjust_text
from radio_module import *
from rotation_script import *
from copy import deepcopy

rng = np.random.default_rng()

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

L_EU_1 = L_EU[:4]
L_EU_2 = L_EU[4:]
Freq_1 = Freq[:4]
Freq_2 = Freq[4:]

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

# ----------------------------
# MWA

d_mwa = {"Frequencies": [[72.30, 103.04], [103.04, 133.76], [138.88, 169.60], [169.60, 200.32], [200.32, 231.04]],
         "RMS Noise": [np.array([24, 24]), np.array([14, 14]), np.array([7, 7]), np.array([5, 5]), np.array([5, 5])]}  # mJy, 2min integration time, 40kHz bandwidth

MWA = pd.DataFrame(data=d_mwa)
integration_time = 8 * 60  # minutes
MWA["RMS Noise"] *= np.sqrt((2 / integration_time)) * 10**(-3) * 5

# Retrieve Data
filename = "NASA2903.csv"
df = pd.read_csv(filename, comment="#")
windfile = "wind_info" + filename[4:-4] + ".txt"

headers_to_find = ["pl_name", "pl_orbper", "pl_orbsmax", "pl_radj", "pl_bmassj", "pl_bmassprov", "pl_dens", "st_rad", "st_mass", "st_age", "sy_dist"]
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
IsBurst = 0

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

y_mag_minerr = []
y_mag_maxerr = []
y_kin_minerr = []
y_kin_maxerr = []
y_both_maxerr = []
y_both_minerr = []
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

bins = [np.logspace(-1, 5, 13), np.logspace(-2.5, 1.5, 13), np.logspace(-3, 1.5, 9),
        np.logspace(-1, 1, 9), np.logspace(-2, 1.5, 15), np.logspace(0.5, 3.5, 10)]

ax1.hist(orbs, bins=bins[0], edgecolor="black", color="xkcd:ocean")
ax2.hist(smas, bins=bins[1], edgecolor="black", color="xkcd:ocean")
ax3.hist(ms,  bins=bins[2], edgecolor="black", color="xkcd:ocean")
ax4.hist(Ms,  bins=bins[3], edgecolor="black", color="xkcd:ocean")
ax5.hist(ts,  bins=bins[4], edgecolor="black", color="xkcd:ocean")
ax6.hist(ds, bins=bins[5], edgecolor="black", color="xkcd:ocean")

hist_noir(ax1)
hist_noir(ax2)
hist_noir(ax3)
hist_noir(ax4)
hist_noir(ax5)
hist_noir(ax6)

xlabels = ["Orbital Period (Days)", "Semi-major Axis (AU)", f"Planet Mass ($M_j$)", "Star Mass ($M_\odot$)", "Star Age (Gyr)", "Distance (pc)"]
for i in range(len(axes)):
    axes[i].set_xlabel(xlabels[i])
    axes[i].set_xscale("log")
    axes[i].set_yscale("log")

fig.text(0.02, 0.30, 'Bin Count', va='center', rotation='vertical', fontsize=10)
fig.text(0.02, 0.75, 'Bin Count', va='center', rotation='vertical', fontsize=10)

fig.supylabel(' \n', va='center', rotation='vertical', fontsize=11)
# fig.suptitle('Distribution of Initial Parameters for the Exoplanet Sample', fontsize=13)

plt.show()

plt.rcParams['font.size'] = 9

detectables_mag = []
detectables_kin = []
detectables_both = []

# Calculate Frequencies and intensities
for i, j in df.iterrows():  # The loop that reads exoplanets from NASA file

    name = j[pl_name]

    T_i = j[pl_orbper]
    T_s = (j[pl_orbper+1] - j[pl_orbper+2]) / 2
    if np.isnan(T_s):
        T_s = T_i / 5

    a_i = j[pl_orbsmax]
    a_s = (j[pl_orbsmax+1] - j[pl_orbsmax+2]) / 2
    if np.isnan(a_s):
        a_s = a_i / 5

    R_i = j[radius]
    R_s = (j[radius+1] - j[radius+2]) / 2
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
    p_s = (j[dens+1] - j[dens+2]) / 2
    if np.isnan(p_s):
        p_s = p_i / 5

    M_s_i = j[st_mass]
    M_s_s = (j[st_mass+1] - j[st_mass+2]) / 2
    if np.isnan(M_s_s):
        M_s_s = M_s_i / 5

    t_i = j[st_age]
    t_s = (j[st_age+1] - j[st_age+2]) / 2
    if np.isnan(t_s):
        t_s = t_i / 5

    Rs_i = j[st_rad]
    Rs_s = (j[st_rad+1] - j[st_rad+2]) / 2
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

    for k in range(1000):  # The loop for Monte Carlo iterations
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
        # d *= 3.26156

        p_c = density(p)
        w_p = rotation(T, a, L, M, R)

        r_c = convective_radius(M, p_c, R)
        mu = magnetic_moment(p_c, w_p, r_c, sigma)
        B = magnetic_field(mu, R)

        if B == 0:
            continue

        D = d * 9.46 * 10 ** 15  # conversion to meters

        # br = B_

        B_perp = imf_perp_complete(M_s, a, Rs, t, v_wind, b0_exponent)
        v_k = keplerian(M_s, a)
        veff = v_eff(v_wind, v_k)
        B_star = imf_complete(M_s, a, v_wind, t, Rs, b0_exponent)
        n = number_density(Mdot, veff, a)
        R_m = Rm(B, R, n, T_wind, v_wind, B_star)

        nu = max_freq(B)
        assert nu > 0, f"Maximum emission frequency must be positive, instead got {nu=}."

        nu /= 10 ** 6

        flag = True

        if freq_condition(nu, n):
            flag = False
            freqs.append(nu)
            nu *= 10**6

            I = complete(B, a, M_s, Mdot, D)
            assert I > 0, f"Radio brightness must be positive, instead got {I=}."

            P_in_mag = P_input_mag(B_perp, veff*10**3, R_m*7*10**8, n)
            P_in_kin = P_input_kin(B_perp, veff*10**3, R_m*7*10**8, n)
            P_rad_mag = radio_power(P_in_mag, 0, nu, D)
            P_rad_kin = radio_power(0, P_in_kin, nu, D)
            P_rad_both = radio_power(P_in_mag, P_in_kin, nu, D, both=True)

            if IsBurst:
                I = I * (10 ** 1.53)
                P_rad_mag *= 10**1.53
                P_rad_kin *= 10**1.53
                P_rad_both *= 10**1.53

            intens_mag.append(P_rad_mag)
            intens_kin.append(P_rad_kin)
            intens_both.append(P_rad_both)

    if flag:
        continue

    # print(f"{B_perp=}, {B_star=}, {n=}, {R_m=}. {Mdot=}, {P_in}, {P_rad}")
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

    obs_mag = ""
    obs_kin = ""
    obs_both = ""

    EXO = Exoplanet(name, a, R, M, p, B, M_s, Mdot, d, freq=nu, intensity_mag=I_mag, intensity_kin=I_kin, intensity_both=I_both)
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
    if 30 <= nu <= 75:
        for m in range(4):
            if Freq_1[m] <= nu <= Freq_1[m + 1]:
                if I_mag >= (L_EU_1[m + 1] - L_EU_1[m]) / (Freq_1[m + 1] - Freq_1[m]) * (nu - Freq_1[m]) + L_EU_1[m]:
                    obs_mag = str(EXO.name)
                    print(obs_mag)
                    detectables_mag.append(EXO)
                    observable_mag = True
                    observable_flag = True
                if I_kin >= (L_EU_1[m + 1] - L_EU_1[m]) / (Freq_1[m + 1] - Freq_1[m]) * (nu - Freq_1[m]) + L_EU_1[m]:
                    obs_kin = str(EXO.name)
                    print(obs_kin)
                    detectables_kin.append(EXO)
                    observable_flag = True
                    observable_kin = True
                if I_both >= (L_EU_1[m + 1] - L_EU_1[m]) / (Freq_1[m + 1] - Freq_1[m]) * (nu - Freq_1[m]) + L_EU_1[m]:
                    obs_both = str(EXO.name)
                    print(obs_both)
                    detectables_both.append(EXO)
                    observable_both = True
                    observable_flag = True

    if 120 <= nu <= 180:
        for m in range(4, 7):
            if Freq_2[m] <= nu <= Freq_2[m+1]:
                if I_mag >= (L_EU_2[m + 1] - L_EU_2[m]) / (Freq_2[m + 1] - Freq_2[m]) * (nu - Freq_2[m]) + L_EU_2[m]:
                    obs_mag = str(EXO.name)
                    print(obs_mag)
                    detectables_mag.append(EXO)
                    observable_mag = True
                if I_kin >= (L_EU_2[m + 1] - L_EU_2[m]) / (Freq_2[m + 1] - Freq_2[m]) * (nu - Freq_2[m]) + L_EU_2[m]:
                    obs_kin = str(EXO.name)
                    print(obs_kin)
                    detectables_kin.append(EXO)
                    observable_kin = True
                if I_both >= (L_EU_2[m + 1] - L_EU_2[m]) / (Freq_2[m + 1] - Freq_2[m]) * (nu - Freq_2[m]) + L_EU_2[m]:
                    obs_both = str(EXO.name)
                    print(obs_both)
                    detectables_both.append(EXO)
                    observable_both = True

    if 72.30 <= nu <= 231.04:
        for m in range(5):
            if not observable_mag:
                if MWA["Frequencies"][m][0] <= nu <= MWA["Frequencies"][m][1] and I_mag > MWA["RMS Noise"][m][0]:
                    obs_mag = str(EXO.name)
                    print(obs_mag)
                    detectables_mag.append(EXO)
                    observable_mag = True

            if not observable_kin:
                if MWA["Frequencies"][m][0] <= nu <= MWA["Frequencies"][m][1] and I_kin > MWA["RMS Noise"][m][0]:
                    obs_kin = str(EXO.name)
                    print(obs_kin)
                    detectables_kin.append(EXO)
                    observable_kin = True

            if not observable_both:
                if MWA["Frequencies"][m][0] <= nu <= MWA["Frequencies"][m][1] and I_both > MWA["RMS Noise"][m][0]:
                    obs_both = str(EXO.name)
                    print(obs_both)
                    detectables_both.append(EXO)
                    observable_both = True
                    break


    # if not observable_mag and not observable_kin and not observable_both:
    if 120 <= nu < 850 or 1050 < nu <= 1450:
        for m in range(4):

            if not observable_mag:
                if uGMRT["Frequencies"][m][0] <= nu <= uGMRT["Frequencies"][m][1] and I_mag > uGMRT["RMS Noise"][m][0]:
                    obs_mag = str(EXO.name)
                    print(obs_mag)
                    detectables_mag.append(EXO)
                    observable_mag = True

            if not observable_kin:
                if uGMRT["Frequencies"][m][0] <= nu <= uGMRT["Frequencies"][m][1] and I_kin > uGMRT["RMS Noise"][m][0]:
                    obs_kin = str(EXO.name)
                    print(obs_kin)
                    detectables_kin.append(EXO)
                    observable_kin = True

            if not observable_both:
                if uGMRT["Frequencies"][m][0] <= nu <= uGMRT["Frequencies"][m][1] and I_both > uGMRT["RMS Noise"][m][0]:
                    obs_both = str(EXO.name)
                    print(obs_both)
                    detectables_both.append(EXO)
                    observable_both = True
                    break

    if EXO.name == "tau Boo b":  # Special interest
        obs = EXO.name

    labels_mag.append(obs_mag)
    labels_kin.append(obs_kin)
    labels_both.append(obs_both)

    labels = [labels_mag, labels_kin, labels_both]

    selected = "tau Boo b"
    if EXO.name == selected:
        plt.subplot(1, 2, 1)
        plt.hist(intens_mag, edgecolor="black")
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

# arr = np.array(labels)
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

df = pd.DataFrame({"x": frequencies,
                    "y_mag": intensities_mag,
                    "y_kin": intensities_kin,
                    "y_both": intensities_both,
                    "d": distances,
                    "s": semis,
                    "l_mag": labels_mag,
                    "l_kin": labels_kin,
                    "l_both": labels_both})

det_data = [[exo.name, exo.freq, exo.intensity_both] for exo in detectables_both]
df_det = pd.DataFrame(det_data[1:], columns=["Name", "Freq", "Flux"])


def scatter_plot(df, which, y_err, x_err, det, zoom=False, save=False, fix_lim=False):

    rc = {"font.family": "sans-serif", "font.weight": "light", "font.variant": "small-caps", "font.size": 10}
    plt.rcParams.update(rc)

    plt.rcParams['figure.figsize'] = [10, 5]

    y_mag_err, y_kin_err, y_both_err = y_err[0], y_err[1], y_err[2]

    df["labels"] = df.apply(lambda row: str(row['l_mag']) + row['l_kin'] + row["l_both"], axis=1)

    cond = df["labels"][df["labels"] == ""].index
    cond_mag = df["l_mag"][df["l_mag"] == ""].index
    cond_kin = df["l_kin"][df["l_kin"] == ""].index
    cond_both = df["l_both"][df["l_both"] == ""].index

    y_mag_err[0][cond_mag] = 0
    y_mag_err[1][cond_mag] = 0
    y_kin_err[0][cond_kin] = 0
    y_kin_err[1][cond_kin] = 0
    y_both_err[0][cond_both] = 0
    y_both_err[1][cond_both] = 0

    x_mag_err = deepcopy(x_err)
    x_kin_err = deepcopy(x_err)
    x_both_err = deepcopy(x_err)
    x_mag_err[0][cond_mag] = 0
    x_mag_err[1][cond_mag] = 0
    x_kin_err[0][cond_kin] = 0
    x_kin_err[1][cond_kin] = 0
    x_both_err[0][cond_both] = 0
    x_both_err[1][cond_both] = 0

    # x_err[0][cond] = 0
    # x_err[1][cond] = 0

    fig0, ax0 = plt.subplots()

    if which == "mag":
        y_err = [y_mag_minerr, y_mag_maxerr]
        x_err_new = x_mag_err
        df["y"] = df["y_mag"]

    elif which == "kin":
        y_err = [y_kin_minerr, y_kin_maxerr]
        x_err_new = x_kin_err
        df["y"] = df["y_kin"]

    else:
        y_err = [y_both_minerr, y_both_maxerr]
        x_err_new = x_both_err
        df["y"] = df["y_both"]

    if zoom:
        size = 20
    else:
        size = df.d

    im = ax0.scatter(df.x, df.y, c=df.s, s=size, cmap="magma_r")
    errorbar = ax0.errorbar(df.x, df.y,
                     yerr=y_err,
                     xerr=x_err_new,
                     fmt="None",
                     ecolor="black",
                     elinewidth=0.5)
    ax0.plot(Freq_1, L_EU_1, linestyle="dashed", color="red", linewidth=0.5)
    ax0.plot(Freq_2, L_EU_2, linestyle="dashed", color="purple", linewidth=0.5)
    for i in range(4):
        x = uGMRT["Frequencies"][i]
        y = uGMRT["RMS Noise"][i]
        plt.plot(x, y, "b--", linewidth=0.5)
        if i == 0:
            ax0.fill_between(x, y, 10 ** 6, color="blue", alpha=0.1, label="uGMRT")
        else:
            ax0.fill_between(x, y, 10 ** 6, color="blue", alpha=0.1)

    for i in range(5):
        x = MWA["Frequencies"][i]
        y = MWA["RMS Noise"][i]
        plt.plot(x, y, "k--", linewidth=0.5)
        if i == 0:
            ax0.fill_between(x, y, 10 ** 6, color="grey", alpha=0.1, label="MWA")
        else:
            ax0.fill_between(x, y, 10 ** 6, color="grey", alpha=0.1)

    ax0.fill_between(Freq_1, L_EU_1, 10**6, color="red", alpha=0.1, label="LOFAR LBA")
    ax0.fill_between(Freq_2, L_EU_2, 10**6, color="purple", alpha=0.1, label="LOFAR HBA")
    plt.colorbar(im, ax=ax0, label="Distance to Host Star ($\log_{10}{\mathrm{(AU)}}$)", aspect=25, extend="both")
    ax0.axvline(x=10, color="black", linestyle="dashed")
    ax0.set_xscale("log")
    ax0.set_yscale("log")
    ax0.axvspan(0, 10, alpha=0.2, color="teal")

    if which == "mag":
        lab = df["l_mag"]
        tit = "\n(Magnetic Energy)"
    elif which == "kin":
        lab = df["l_kin"]
        tit = "\n(Kinetic Energy)"
    else:
        lab = df["l_both"]
        tit = ""

    df["x_errmin"], df["x_errmax"] = x_err_new
    df["y_errmin"], df["y_errmax"] = y_err

    if zoom:
        errorbar.remove()
        df1 = df[lab != ""]
        yerr = [df1["y_errmin"], df1["y_errmax"]]
        xerr = [df1["x_errmin"], df1["x_errmax"]]
        ax0.errorbar(df1.x, df1.y, yerr=yerr, xerr=xerr, fmt="None", ecolor="black", elinewidth=1, capsize=2)
        ax0.set_xlim(left=min(det["Freq"]) * 0.5, right=max(det["Freq"]) * 2)
        ax0.set_ylim(bottom=min(det["Flux"]) * 0.5, top=max(det["Flux"]) * 2)
        texts = [plt.text(df.x[i], df.y[i], lab[i], ha='center', va='center', fontsize=8) for i in range(len(lab)) if lab[i] != ""]
        adjust_text(texts, arrowprops=dict(arrowstyle="-", color="r", lw=0.5))
        plt.legend()

    else:
        ax0.set_xlim(left=0.05)
        if fix_lim:
            ax0.set_ylim(bottom=1e-10, top=1)
        else:
            ax0.set_ylim(bottom=min(df.y) * 0.1, top=max(df.y) * 10)
        # for i, txt in enumerate(labels):
        #     if txt:
        #         ax0.annotate(txt, xy=(df1.x[i], df1.y[i]), xytext=(2, 2), textcoords="offset pixels", fontsize=7)
        plt.legend(loc="upper left")

    # if IsBurst:
    #     ax0.set_title(f"Frequency and Intensity of Burst CMI Emissions of the Exoplanet Sample{tit}")
    # else:
    #     ax0.set_title(f"Frequency and Intensity of Quiescent CMI Emissions of the Exoplanet Sample{tit}")

    ax0.set_xlabel("Emission Frequency (MHz)")
    ax0.set_ylabel("Radio Brightness (Jy)")
    retro_noir(ax0)
    fig0.tight_layout()

    if save:
        if zoom:
            plt.savefig("zoom.pdf")
        else:
            plt.savefig("scatter.pdf")


def outcome_dist_hists(intensities, which, magnetic_fields, save=False):
    if which == "mag":
        intensities = intensities[0]
    elif which == "kin":
        intensities = intensities[1]
    else:
        intensities = intensities[2]

    bin1_lower = math.floor(math.log10(min(intensities)))
    bin1_higher = math.floor(math.log10(max(intensities))) + 1
    n = (bin1_higher - bin1_lower) + 1
    TheBins1 = np.logspace(bin1_lower, bin1_higher, n)

    plt.rcParams['figure.figsize'] = [6, 4]
    rc = {"font.size": 10}
    plt.rcParams.update(rc)

    fig1, axs = plt.subplots(1, 2, sharey="row", figsize=[10, 5])

    ax1, ax2 = axs[0], axs[1]

    ax1.hist(intensities, bins=TheBins1, edgecolor="black", color="xkcd:sea")
    if IsBurst:
        ax1.set_xlabel("Intensity of Burst Emission (Jy)")
        # ax1.set_title("Histogram of Burst Emission Intensities")
    else:
        ax1.set_xlabel("Intensity of Quiescent Emission (Jy)")
        # ax1.set_title("Histogram of Quiescent Emission Intensities")

    ax1.set_xscale("log")
    ax1.set_yscale("log")

    bin2_lower = math.floor(math.log10(min(magnetic_fields)))
    bin2_higher = math.floor(math.log10(max(magnetic_fields)))
    n = (bin2_higher - bin2_lower) * 2 + 1
    TheBins2 = np.logspace(bin2_lower, bin2_higher, n)

    ax2.hist(magnetic_fields, bins=TheBins2, edgecolor="black", color="xkcd:sea")
    ax2.set_xlabel("Magnetic Field Strength at the Surface (Gauss)")
    # ax2.set_title("Histogram of the Magnetic Field Strengths")

    ax2.set_xscale("log")
    hist_noir(ax1)
    hist_noir(ax2)
    fig1.supylabel("Number of Exoplanets")

    if save:
        plt.savefig("hist.pdf")


scatter_plot(df, "both", y_err, x_err, df_det)
outcome_dist_hists(intensities, "both", magnetic_fields)

plt.show()

for i in range(len(intensities)):
    l1 = [frequencies, intensities[i]]
    l1 = [arr.tolist() for arr in l1]
    l2 = [names]
    l2.extend(l1)

    if i == 0:
        file_names = ["names_mag.txt", "freq_mag.txt", "intens_mag.txt", "detectables_mag.csv"]
        detectables = detectables_mag
        detectable_data = pd.DataFrame([[exo.name, exo.freq, exo.intensity_mag] for exo in detectables], columns=["Name", "Freq", "Flux"])

    elif i == 1:
        file_names = ["names_kin.txt", "freq_kin.txt", "intens_kin.txt", "detectables_kin.csv"]
        detectables = detectables_kin
        detectable_data = pd.DataFrame([[exo.name, exo.freq, exo.intensity_kin] for exo in detectables], columns=["Name", "Freq", "Flux"])

    else:
        file_names = ["names_both.txt", "freq_both.txt", "intens_both.txt", "detectables_both.csv"]
        detectables = detectables_both
        detectable_data = pd.DataFrame([[exo.name, exo.freq, exo.intensity_both] for exo in detectables], columns=["Name", "Freq", "Flux"])

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

# det = np.genfromtxt("detectables_both.txt", usecols=0, dtype=str, delimiter="  ", skip_header=1, skip_footer=1)
# detect_data = df[df["pl_name"].isin(det)]
# enough_data = detect_data[["pl_name", "pl_bmassj", "pl_radj", "pl_orbsmax", "sy_dist", "st_age"]]
# sorted = enough_data.sort_values(by="pl_name")


import numpy as np
import pandas as pd
import math
from tabulate import tabulate
from adjustText import adjust_text
from radio_module import *
from rotation_script import *
from copy import deepcopy
import smplotlib
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D



def attempt_hershley():
    from fontTools.ttLib import TTFont
    from matplotlib.font_manager import FontProperties
    from matplotlib import font_manager
    from matplotlib import mathtext

    font_path = '/Users/asafkaya/Documents/My Stuff/Programming/PythonFiles/Hershey_font_TTF-main/ttf/AVHersheyComplexHeavy.ttf'
    font_manager.fontManager.addfont(font_path)
    font_prop = font_manager.FontProperties(fname=font_path)
    custom_font_name = font_prop.get_name()
    plt.rcParams['font.family'] = custom_font_name
    plt.rcParams['font.weight'] = "heavy"

    def replace_minus_with_hyphen(s):
        return s.replace('\u2212', '\u002D')

    # Override the default text rendering to replace minus signs globally

    class CustomScalarFormatter(mpl.ticker.ScalarFormatter):
        def format_data(self, value):
            formatted_value = super().format_data(value)
            return replace_minus_with_hyphen(formatted_value)

    # Apply the custom formatter to tick labels
    def apply_custom_formatter(ax):
        ax.xaxis.set_major_formatter(CustomScalarFormatter())
        ax.yaxis.set_major_formatter(CustomScalarFormatter())

    # Custom math text rendering
    mpl.rcParams['text.usetex'] = True
    mpl.rcParams['text.latex.preamble'] = r'\usepackage{textcomp}'


plt.rcParams['figure.figsize'] = [10, 5]

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
# NenuFAR
NenuNoise = np.array([130.0, 9.0])  # MJy (10 MHz, 1h, website)
NenuFreq = np.array([15.0, 85.0])  # MHz (website)
NenuNoise *= np.sqrt(1 / 8) * 10**(-3) * 5  # 5 sigma sensitivity in Jy

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
plt.rcParams['font.size'] = 10

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

# hist_noir(ax1)
# hist_noir(ax2)
# hist_noir(ax3)
# hist_noir(ax4)
# hist_noir(ax5)
# hist_noir(ax6)

xlabels = ["Orbital Period (Days)", "Semi-major Axis (AU)", f"Planet Mass ($M_j$)", "Star Mass ($M_\odot$)", "Star Age (Gyr)", "Distance (pc)"]
for i in range(len(axes)):
    axes[i].set_xlabel(xlabels[i])
    axes[i].set_xscale("log")
    axes[i].set_yscale("log")

fig.text(0.04, 0.30, 'Bin Count', va='center', rotation='vertical', fontsize=12)
fig.text(0.04, 0.75, 'Bin Count', va='center', rotation='vertical', fontsize=12)

fig.tight_layout()

fig.supylabel(' \n', va='center', rotation='vertical', fontsize=11)
# fig.suptitle('Distribution of Initial Parameters for the Exoplanet Sample', fontsize=13)


plt.show()

# plt.rcParams['font.size'] = 9

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
    n_ps = []

    # high_var = ["AU Mic c", "V1298 Tau b"]
    high_var = ["AU Mic c", "V1298 Tau d", "V1298 Tau b", "V1298 Tau e", "V1298 Tau c"]
    # high_var = []

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

        n_p = 8.98 * np.sqrt(n) * 10**(-3)

        nu = max_freq(B)
        assert nu > 0, f"Maximum emission frequency must be positive, instead got {nu=}."

        nu /= 10 ** 6
        n_p /= 10**6

        flag = True

        if freq_condition(nu, n) or name in high_var:
            flag = False
            freqs.append(nu)
            n_ps.append(n_p)
            nu *= 10**6

            I = complete(B, a, M_s, Mdot, D)
            assert I > 0, f"Radio brightness must be positive, instead got {I=}."

            P_in_mag = P_input_mag(B_perp, veff*10**3, R_m*7*10**8, n)
            P_in_kin = P_input_kin(B_perp, veff*10**3, R_m*7*10**8, n)
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

    # print(f"{B_perp=}, {B_star=}, {n=}, {R_m=}. {Mdot=}, {P_in}, {P_rad}")
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

    x_low.append((np.percentile(freqs, 16))/nu)
    x_up.append((np.percentile(freqs, 84))/nu)

    y_low.append((np.percentile(intens_both, 16))/I_both)
    y_up.append((np.percentile(intens_both, 84))/I_both)

    x_err_avg = (-np.log10((np.percentile(freqs, 16))/nu) + np.log10((np.percentile(freqs, 84))/nu)) / 2
    y_err_avg = (-np.log10((np.percentile(intens_both, 16))/I_both) + np.log10((np.percentile(intens_both, 84))/I_both)) / 2

    obs_mag = ""
    obs_kin = ""
    obs_both = ""
    out = ""

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
    observable_flag = False
    # if 30 <= nu <= 75:
    # bw = 4.66  # LOFAR bandwidth of 4.66 MHz

    if x_err_avg < 0.85:

        for m in range(3):
            x = freqs[(Freq_1[m] <= freqs) & (freqs <= Freq_1[m + 1])]
            frac = len(x) / len(freqs)

            if frac > 0.1:
                # if I_mag >= (L_EU_1[m + 1] - L_EU_1[m]) / (Freq_1[m + 1] - Freq_1[m]) * (nu - Freq_1[m]) + L_EU_1[m]:
                if I_mag >= (L_EU_1[m] + L_EU_1[m+1]) / 2:
                    obs_mag = str(EXO.name)
                    print(obs_mag, m, I_mag)
                    observable_mag = True
                # if I_kin >= (L_EU_1[m + 1] - L_EU_1[m]) / (Freq_1[m + 1] - Freq_1[m]) * (nu - Freq_1[m]) + L_EU_1[m]:
                if I_kin >= (L_EU_1[m] + L_EU_1[m+1]) / 2:
                    obs_kin = str(EXO.name)
                    print(obs_kin, m, I_kin)
                    observable_kin = True
                # if I_both >= (L_EU_1[m + 1] - L_EU_1[m]) / (Freq_1[m + 1] - Freq_1[m]) * (nu - Freq_1[m]) + L_EU_1[m]:
                if I_both >= (L_EU_1[m] + L_EU_1[m+1]) / 2:
                    obs_both = str(EXO.name)
                    print(obs_both, m, I_both)
                    observable_both = True
                    lofar_obs.append(obs_both)
                    if not Freq_1[0] <= nu <= Freq_1[3]:
                        out = str(EXO.name)
                        print(f"{out} Outlier for LOFAR LBA {m=}")
                        outliers.add(obs_both)
                    else:
                        insiders.add(obs_both)
                observable_flag = observable_both or observable_mag or observable_kin

        # if 120 <= nu <= 180:
        # bw = 4.66  # LOFAR bandwidth of 4.66 MHz

        for m in range(4, 6):
            x = freqs[(Freq_2[m] <= freqs) & (freqs <= Freq_2[m + 1])]
            frac = len(x) / len(freqs)

            if frac > 0.1:
                # if I_mag >= (L_EU_2[m + 1] - L_EU_2[m]) / (Freq_2[m + 1] - Freq_2[m]) * (nu - Freq_2[m]) + L_EU_2[m]:
                if I_mag >= (L_EU_2[m] + L_EU_2[m+1]) / 2:
                    obs_mag = str(EXO.name)
                    print(obs_mag)
                    observable_mag = True

                # if I_kin >= (L_EU_2[m + 1] - L_EU_2[m]) / (Freq_2[m + 1] - Freq_2[m]) * (nu - Freq_2[m]) + L_EU_2[m]:
                if I_kin >= (L_EU_2[m] + L_EU_2[m+1]) / 2:
                    obs_kin = str(EXO.name)
                    print(obs_kin)
                    observable_kin = True
                # if I_both >= (L_EU_2[m + 1] - L_EU_2[m]) / (Freq_2[m + 1] - Freq_2[m]) * (nu - Freq_2[m]) + L_EU_2[m]:

                if I_both >= (L_EU_2[m] + L_EU_2[m+1]) / 2:
                    obs_both = str(EXO.name)
                    print(obs_both)
                    observable_both = True
                    lofar_obs.append(obs_both)
                    if not Freq_2[4] <= nu <= Freq_2[6]:
                        out = str(EXO.name)
                        print(f"{out} Outlier for LOFAR HBA {m=} ")
                        outliers.add(obs_both)
                    else:
                        insiders.add(obs_both)
                observable_flag = observable_both or observable_mag or observable_kin

        # NenuFAR
        x = freqs[(NenuFreq[0] <= freqs) & (freqs <= NenuFreq[1])]
        frac = len(x) / len(freqs)

        # if frac > 0.1:

        if NenuFreq[0] <= nu <= NenuFreq[1]:
            if I_mag >= NenuNoise[0] + (NenuNoise[1] - NenuNoise[0]) / (NenuFreq[1] - NenuFreq[0]) * (nu - NenuFreq[0]):
                obs_mag = str(EXO.name)
                print(obs_mag)
                observable_mag = True

            if I_kin >= NenuNoise[0] + (NenuNoise[1] - NenuNoise[0]) / (NenuFreq[1] - NenuFreq[0]) * (nu - NenuFreq[0]):
                obs_kin = str(EXO.name)
                print(obs_kin)
                observable_kin = True

            if I_both >= NenuNoise[0] + (NenuNoise[1] - NenuNoise[0]) / (NenuFreq[1] - NenuFreq[0]) * (nu - NenuFreq[0]):
                obs_both = str(EXO.name)
                print(obs_both)
                observable_both = True
                nenufar_obs.append(obs_both)
                insiders.add(obs_both)

        elif frac > 0.1 and I_both > (NenuNoise[0] + NenuNoise[1]) / 2:
            obs_both = str(EXO.name)
            print(obs_both)
            observable_both = True
            nenufar_obs.append(obs_both)
            out = str(EXO.name)
            print(f"{out} Outlier for NenuFAR")
            outliers.add(obs_both)
        # observable_flag = observable_both or observable_mag or observable_kin


        # if 72.30 <= nu <= 231.04:
        for m in range(5):
            x = freqs[(MWA["Frequencies"][m][0] <= freqs) & (freqs <= MWA["Frequencies"][m][1])]
            frac = len(x) / len(freqs)

            if frac > 0.1 and I_mag > MWA["RMS Noise"][m][0]:
                obs_mag = str(EXO.name)
                print(obs_mag)
                observable_mag = True

            if frac > 0.1 and I_kin > MWA["RMS Noise"][m][0]:
                obs_kin = str(EXO.name)
                print(obs_kin)
                observable_kin = True

            if frac > 0.1 and I_both > MWA["RMS Noise"][m][0]:
                obs_both = str(EXO.name)
                print(obs_both)
                observable_both = True
                mwa_obs.append(obs_both)
                if not MWA["Frequencies"][0][0] <= nu <= MWA["Frequencies"][4][1]:
                    out = str(EXO.name)
                    print(f"{out} Outlier for MWA {m=} ")
                    outliers.add(obs_both)
                else:
                    insiders.add(obs_both)
                break

        for m in range(4):

            x = freqs[(uGMRT["Frequencies"][m][0] <= freqs) & (freqs <= uGMRT["Frequencies"][m][1])]
            frac = len(x) / len(freqs)

            if frac > 0.1 and I_mag > uGMRT["RMS Noise"][m][0]:
                obs_mag = str(EXO.name)
                print(obs_mag)
                observable_mag = True

            if frac > 0.1 and I_kin > uGMRT["RMS Noise"][m][0]:
                obs_kin = str(EXO.name)
                print(obs_kin)
                observable_kin = True

            if frac > 0.1 and I_both > uGMRT["RMS Noise"][m][0]:
                obs_both = str(EXO.name)
                print(obs_both)
                observable_both = True
                ugmrt_obs.append(obs_both)
                if not uGMRT["Frequencies"][0][0] <= nu <= uGMRT["Frequencies"][3][1]:
                    out = str(EXO.name)
                    print(f"{out} Outlier for uGMRT {m=}")
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
        ax1.axvline(np.log10(np.percentile(intens_both, 16)), color="xkcd:deep green", linestyle="--", label="16th & 84th\npercentiles")
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
        fig.legend(lines, labels, loc="upper right", bbox_to_anchor=(1, 1), ncol=int(len(lines_labels)/2), frameon=True, shadow=True)
        plt.savefig(f"{EXO.name}.pdf")
        # plt.show()
        # plt.close()

# arr = np.array(labels)
real_outliers = set(outliers) - set(insiders)
lofar_obs, nenufar_obs, mwa_obs, ugmrt_obs = np.array(lofar_obs), np.array(nenufar_obs), np.array(mwa_obs), np.array(ugmrt_obs)
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

det_data = [[exo.name, exo.freq, exo.intensity_both] for exo in detectables_both]
df_det = pd.DataFrame(det_data[0:], columns=["Name", "Freq", "Flux"])


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
    # rc = {"font": font}
    rc = {"font.size": 12}
    plt.rcParams.update(rc)

    fig1, axs = plt.subplots(1, 2, sharey="row", figsize=[10, 5])

    ax1, ax2 = axs[0], axs[1]

    ax1.hist(intensities, bins=TheBins1, edgecolor="black", color="xkcd:sea")
    if IsBurst:
        ax1.set_xlabel("Flux Density of Burst Emission (Jy)")
        # ax1.set_title("Histogram of Burst Emission Intensities")
    else:
        ax1.set_xlabel("Flux Density of Quiescent Emission (Jy)")
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
    # hist_noir(ax1)
    # hist_noir(ax2)
    fig1.supylabel("Number of Exoplanets")

    if save:
        plt.savefig("hist.pdf")


def is_within_limits(x, y, xlim, ylim):
    return xlim[0] <= x <= xlim[1] and ylim[0] <= y <= ylim[1]


# The following function is a truly badly-written one. I have stopped caring for its readability at this point. Sorry about this. At least it gets the job done.
def scatter_plot(df1, which, y_err, x_err, det, avg_err, zoom=False, save=False, fix_lim=False, strict=False, others=0):

    df = df1.copy()

    plt.rcParams['figure.figsize'] = [10, 5]
    plt.rcParams['font.size'] = 12

    y_mag_err, y_kin_err, y_both_err, y_strict_err = y_err[0], y_err[1], y_err[2], y_err[2]

    df["labels"] = df.apply(lambda row: str(row['l_mag']) + row['l_kin'] + row["l_both"], axis=1)

    cond = df["labels"][df["labels"] == ""].index
    cond_mag = df["l_mag"][df["l_mag"] == ""].index
    cond_kin = df["l_kin"][df["l_kin"] == ""].index
    cond_both = df["l_both"][df["l_both"] == ""].index
    cond_strict = df["l_both"][(df["l_both"] == "") & (df["l_both"].isin(real_outliers))].index


    y_mag_err[0][cond_mag] = 0
    y_mag_err[1][cond_mag] = 0
    y_kin_err[0][cond_kin] = 0
    y_kin_err[1][cond_kin] = 0
    y_both_err[0][cond_both] = 0
    y_both_err[1][cond_both] = 0
    y_strict_err[0][cond_strict] = 0
    y_strict_err[1][cond_strict] = 0


    x_mag_err = deepcopy(x_err)
    x_kin_err = deepcopy(x_err)
    x_both_err = deepcopy(x_err)
    x_strict_err = deepcopy(x_err)
    x_mag_err[0][cond_mag] = 0
    x_mag_err[1][cond_mag] = 0
    x_kin_err[0][cond_kin] = 0
    x_kin_err[1][cond_kin] = 0
    x_both_err[0][cond_both] = 0
    x_both_err[1][cond_both] = 0
    x_strict_err[0][cond_strict] = 0
    x_strict_err[1][cond_strict] = 0

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

    df["xerr_nonzero"] = np.where(df["x"] == 0, np.nan, df["x"])
    df["yerr_nonzero"] = np.where(df["y"] == 0, np.nan, df["y"])

    if which == "mag":
        lab = df["l_mag"]
        tit = "\n(Magnetic Energy)"
    elif which == "kin":
        lab = df["l_kin"]
        tit = "\n(Kinetic Energy)"
    else:
        lab = df["l_both"]
        tit = ""

    lab_strict = df["l_both"].copy()
    lab_strict[df["l_both"].isin(real_outliers)] = ""

    if zoom:
        if strict:
            size = lab_strict.apply(lambda x: others if x == "" else 60)
        else:
            size = lab.apply(lambda x: others if x == "" else 60)

    else:
        size = df.d

    if not zoom:
        smplotlib.set_style(edgecolor='face')
    else:
        smplotlib.set_style(edgecolor='k')

    outliers = df[df['Names'].isin(real_outliers)]
    non_outliers = df[~df['Names'].isin(real_outliers)]
    min_color = df['s'].min()
    max_color = df['s'].max()

    if not zoom:

        scatter_non_outliers = ax0.scatter(non_outliers.xerr_nonzero, non_outliers.yerr_nonzero, s=non_outliers.d,
                                          c=non_outliers.s, cmap='magma_r', marker='o', norm=Normalize(vmin=min_color, vmax=max_color))

        # Plot outliers with upside-down triangles, using size and color mappings
        scatter_outliers = ax0.scatter(outliers.xerr_nonzero, outliers.yerr_nonzero, s=outliers.d,
                                      c=outliers.s, cmap='magma_r', marker='v', norm=Normalize(vmin=min_color, vmax=max_color))

    else:
        scatter_non_outliers = ax0.scatter(non_outliers.xerr_nonzero, non_outliers.yerr_nonzero, s=size[~df['Names'].isin(real_outliers)],
                                           c=non_outliers.s, cmap='magma_r', marker='o', norm=Normalize(vmin=min_color, vmax=max_color))

        # Plot outliers with upside-down triangles, using size and color mappings
        scatter_outliers = ax0.scatter(outliers.xerr_nonzero, outliers.yerr_nonzero, s=size[df['Names'].isin(real_outliers)],
                                       c=outliers.s, cmap='magma_r', marker='v', norm=Normalize(vmin=min_color, vmax=max_color))

    # im = ax0.scatter(df.xerr_nonzero, df.yerr_nonzero, c=df.s, s=size, cmap="magma_r")
    errorbar = ax0.errorbar(df.x, df.y,
                     yerr=y_err,
                     xerr=x_err_new,
                     fmt="None",
                     ecolor="black",
                     elinewidth=0.5,
                     capsize=0)
    # errorbar = ax0.errorbar(df.x, df.y,
    #                  yerr=y_err_clean,
    #                  xerr=x_err_clean,
    #                  fmt="None",
    #                  ecolor="black",
    #                  elinewidth=0.5,
    #                  capsize=0)

    x_lba = []
    y_lba = []
    for i in range(len(Freq_1) - 1):
        x = np.linspace(Freq_1[i], Freq_1[i + 1], 50).tolist()
        y = np.linspace(L_EU_1[i], L_EU_1[i + 1], 50).tolist()
        x_lba.extend(x)
        y_lba.extend(y)

    ax0.plot(x_lba, y_lba, linestyle="-", color="red", linewidth=0.5)
    ax0.fill_between(x_lba, y_lba, 10 ** 6, color="red", alpha=0.1, label="LOFAR LBA")

    x_hba = []
    y_hba = []
    for i in range(len(Freq_1), (len(Freq_1) + len(Freq_2) - 1)):
        x = np.linspace(Freq_2[i], Freq_2[i + 1], 50).tolist()
        y = np.linspace(L_EU_2[i], L_EU_2[i + 1], 50).tolist()
        x_hba.extend(x)
        y_hba.extend(y)
    ax0.plot(x_hba, y_hba, linestyle="-", color="purple", linewidth=0.5)
    ax0.fill_between(x_hba, y_hba, 10 ** 6, color="purple", alpha=0.1, label="LOFAR HBA")


    x_nenu = np.linspace(NenuFreq[0], NenuFreq[1], 100)
    y_nenu = np.linspace(NenuNoise[0], NenuNoise[1], 100)
    ax0.plot(x_nenu, y_nenu, "g-", linewidth=0.5)
    ax0.fill_between(x_nenu, y_nenu, 10 ** 6, color="green", alpha=0.1, label="NenuFAR")

    for i in range(4):
        x = uGMRT["Frequencies"][i]
        y = uGMRT["RMS Noise"][i]
        plt.plot(x, y, "b-", linewidth=0.5)
        if i == 0:
            ax0.fill_between(x, y, 10 ** 6, color="blue", alpha=0.1, label="uGMRT")
        else:
            ax0.fill_between(x, y, 10 ** 6, color="blue", alpha=0.1)

    for i in range(5):
        x = MWA["Frequencies"][i]
        y = MWA["RMS Noise"][i]
        plt.plot(x, y, "k-", linewidth=0.5)
        if i == 0:
            ax0.fill_between(x, y, 10 ** 6, color="grey", alpha=0.1, label="MWA")
        else:
            ax0.fill_between(x, y, 10 ** 6, color="grey", alpha=0.1)

    # ax0.fill_between(Freq_1, L_EU_1, 10**6, color="red", alpha=0.1, label="LOFAR LBA")
    # ax0.fill_between(Freq_2, L_EU_2, 10**6, color="purple", alpha=0.1, label="LOFAR HBA")

    norm = Normalize(vmin=df["s"].min(), vmax=df['s'].max())
    sm = ScalarMappable(cmap='magma_r', norm=norm)
    sm.set_array([])  # You can set an array here if needed for specific values

    # cbar = plt.colorbar(im, ax=ax0, label="Distance to Host Star ($\log_{10}{\mathrm{(AU)}}$)", aspect=25, extend="both")
    fig.colorbar(sm, ax=ax0, label="Distance to Host Star ($\log_{10}{\mathrm{(AU)}}$)", aspect=25, extend="both")

    # cbar.ax.tick_params(labelsize=10)
    # cbar.set_label('my label', size='xx-small')
    ax0.axvline(x=10, color="black", linestyle="dashed")
    ax0.set_xscale("log")
    ax0.set_yscale("log")
    ax0.axvspan(0, 10, alpha=0.2, color="teal")

    df["x_errmin"], df["x_errmax"] = x_err_new
    df["y_errmin"], df["y_errmax"] = y_err

    if zoom:
        errorbar.remove()
        ax0.grid("on", alpha=0.2)
        # lab = lab_strict
        if strict:
            df1 = df[lab_strict != ""]
        else:
            df1 = df[lab != ""]
        yerr = [df1["y_errmin"], df1["y_errmax"]]
        xerr = [df1["x_errmin"], df1["x_errmax"]]
        ax0.errorbar(df1.x, df1.y, yerr=yerr, xerr=xerr, fmt="None", ecolor="black", elinewidth=1, capsize=2)
        # ax0.errorbar(df1.x, df1.y, yerr=y_err_inclusive, xerr=x_err_inclusive, fmt="None", ecolor="black", elinewidth=1, capsize=2)
        if fix_lim or strict:
            det = det[~det["Name"].isin(real_outliers)]
        ax0.set_xlim(left=min(det["Freq"]) / 1.25, right=max(det["Freq"]) * 1.25)
        ax0.set_ylim(bottom=min(det["Flux"]) / 1.75, top=max(det["Flux"]) * 1.75)
        xlim = ax0.get_xlim()
        ylim = ax0.get_ylim()
        if strict:
            texts = [plt.text(df.x[i], df.y[i], lab_strict[i], ha='center', va='center', fontsize=8) for i in range(len(lab_strict)) if lab_strict[i] != ""]
        elif fix_lim:
            texts = [plt.text(df.x[i], df.y[i], lab[i], ha='center', va='center', fontsize=8) for i in range(len(lab)) if lab[i] != "" and is_within_limits(df.x[i], df.y[i], xlim, ylim)]
        else:
            texts = [ax0.text(df.x[i], df.y[i], lab[i], ha='center', va='center', fontsize=8) for i in range(len(lab))
                     if lab[i] != ""]
        fig0.legend(fontsize=13, bbox_to_anchor=(0.1, 0.15), loc="lower left", frameon=True)
        line1 = Line2D([0], [0], marker="v", linestyle="None", markerfacecolor="orange", markeredgecolor="black")
        line2 = Line2D([0], [0], marker="o", linestyle="None", markerfacecolor="orange", markeredgecolor="black")
        fig0.legend((line1, line2), ("Outliers", "Insiders"), frameon=True, shadow=True, bbox_to_anchor=(0.8, 0.95), fontsize=12)

        adjust_text(texts, arrowprops=dict(arrowstyle="-", color="k", lw=0.5),
                    force_points=(3, 3), force_text=(2, 2), force_objects=(1.5, 1.5),
                    expand_points=(1.15, 1.15), expand_objects=(1.5, 1.5), expand_align=(1.2, 1.2), precision=20)

    else:
        ax0.set_xlim(left=0.05)
        if fix_lim:
            ax0.set_ylim(bottom=1e-10, top=1)
        else:
            ax0.set_ylim(bottom=min(df.y) * 0.05, top=max(df.y) * 2)
        # for i, txt in enumerate(labels):
        #     if txt:
        #         ax0.annotate(txt, xy=(df1.x[i], df1.y[i]), xytext=(2, 2), textcoords="offset pixels", fontsize=7)
        center_x = 200
        center_y = 3
        arrowprops = dict(arrowstyle='<->,  head_length=0.1', color='k', lw=2)
        ax0.annotate("", xy=(center_x, center_y*avg_err[2]), xytext=(center_x, center_y*avg_err[3]), arrowprops=arrowprops)
        ax0.annotate("", xy=(center_x*avg_err[0], center_y), xytext=(center_x*avg_err[1], center_y), arrowprops=arrowprops)

        ax0.legend(loc="center right", bbox_to_anchor=(1, 0.3), fontsize=11, frameon=True, shadow=True, ncol=1)

        xmin, xmax = ax0.get_xlim()

        # ax0.annotate("Observable From Ground", xy=(np.sqrt(xmax*10), max(df.y)*1.5),
        #     xytext=(np.sqrt(xmax*10), max(df.y)*2), ha='center',  # Horizontal alignment of the text
        #     fontsize=12,  # Size of the text
        #              )

        # ax0.annotate("", xy=(10, max(df.y)*1.5), xytext=(xmax, max(df.y)*1.5), arrowprops=dict(arrowstyle='<-', color='black', linestyle="dotted"))

        ax0.text(np.sqrt(xmax*10), min(df.y)/4, "Observable From Ground", color="green",
            ha='center', va="center",  # Horizontal alignment of the text
            fontsize=12,  # Size of the text
            bbox=dict(facecolor='none', edgecolor='green', boxstyle='square')
                 )

        ax0.text(np.sqrt(xmin*10), min(df.y)/4, "Cannot Penetrate the Ionosphere", color="red",
            ha='center', va="center",  # Horizontal alignment of the text
            fontsize=12,  # Size of the text
            bbox=dict(facecolor='none', edgecolor='red', boxstyle='sawtooth')
                 )

    ax0.set_xlabel("Maximum Emission Frequency (MHz)")
    ax0.set_ylabel("Radio Flux Density (Jy)")
    retro_noir(ax0)
    fig0.tight_layout()

    if save:
        if not zoom:
            plt.savefig("scatter.pdf")
        elif not strict:
            plt.savefig("zoom.pdf")
        elif strict:
            plt.savefig("zoom_inside.pdf")
        elif fix_lim and not strict:
            plt.savefig("zoom_fixed.pdf")


scatter_plot(df1, "both", y_err, x_err, df_det, average_errors, zoom=True)
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
    file_name = f"old_result_tables/{file_names[0]}"
    with open(file_name, 'w') as f:
        f.write(tabulate(table))
        f.close()

    table = sorted(table, key=lambda x: x[1])
    file_name = f"old_result_tables/{file_names[1]}"
    with open(file_name, 'w') as f:
        f.write(tabulate(table))
        f.close()

    table = sorted(table, key=lambda x: x[2], reverse=True)
    file_name = f"old_result_tables/{file_names[2]}"
    with open(file_name, 'w') as f:
        f.write(tabulate(table))
        f.close()

final_df = df1[["Names", "x", "y_both", "x_err_min", "x_err_max", "y_err_min", "y_err_max"]]
final_df.columns = ["Name", "Max. Frequency [MHz]", "Max. Flux Density [Jy]", "Max. Frequency Lower Uncertainty [MHz]", "Max. Frequency Upper Uncertainty [MHz]", "Max. Flux Density Lower Uncertainty [Jy]", "Max. Flux Density Upper Uncertainty [Jy]"]
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

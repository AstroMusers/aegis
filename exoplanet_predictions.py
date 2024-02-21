import pandas as pd
from tabulate import tabulate
from adjustText import adjust_text
from radio_module import *
import rotation_script

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

# Retrieve Data
filename = "NASA3112.csv"
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

bins = [np.logspace(-1, 6, 15), np.logspace(-2.5, 1.5, 13), np.logspace(-3, 1.5, 9),
        np.logspace(-1.5, 1, 11), np.logspace(-2, 1.5, 15), np.logspace(0.5, 3.5, 10)]

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

plt.rcParams['font.size'] = 9

detectables = []

# Calculate Frequencies and intensities
for i, j in df.iterrows():  # The loop that reads exoplanets from NASA file

    name = j[pl_name]

    T_i = j[pl_orbper]
    T_s = (j[pl_orbper+1] - j[pl_orbper+2]) / 2
    if np.isnan(T_s):
        T_s = 0

    a_i = j[pl_orbsmax]
    a_s = (j[pl_orbsmax+1] - j[pl_orbsmax+2]) / 2
    if np.isnan(a_s):
        a_s = 0

    R_i = j[radius]
    R_s = (j[radius+1] - j[radius+2]) / 2
    if np.isnan(R_s):
        R_s = 0

    if j[pl_bmassprov] == "Mass" or j[pl_bmassprov] == "Msin(i)/sin(i)":
        M_i = j[pl_bmassj]
    else:
        M_i = j[pl_bmassj] * 1.15  # Expected Value of the mass based on projected mass
    M_ss = (j[pl_bmassj + 1] - j[pl_bmassj + 2]) / 2
    if np.isnan(M_ss):
        M_ss = 0

    p_i = j[dens]
    p_s = (j[dens+1] - j[dens+2]) / 2
    if np.isnan(p_s):
        p_s = 0

    M_s_i = j[st_mass]
    M_s_s = (j[st_mass+1] - j[st_mass+2]) / 2
    if np.isnan(M_s_s):
        M_s_s = 0

    t_i = j[st_age]
    t_s = (j[st_age+1] - j[st_age+2]) / 2
    if np.isnan(t_s):
        t_s = 0

    Rs_i = j[st_rad]
    Rs_s = (j[st_rad+1] - j[st_rad+2]) / 2
    if np.isnan(Rs_s):
        Rs_s = 0

    T_wind = wind_temperatures[i]
    v_wind = wind_speeds[i]

    d = j[distance]
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

        L = 10 ** rotation_script.kde.resample(1)[0][0]

        highS_Mdot = t ** (-1.23) * 10 ** 3
        lowS_Mdot = t ** (-0.9) * 10 ** 3
        Mdot = 8.1 * t ** (-1.37)  # 10^-14 M_sun / yr = Mdot_sun

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

        B_perp = imf_perp_complete(M_s, a, Rs, t, v_wind)
        v_k = keplerian(M_s, a)
        veff = v_eff(v_wind, v_k)
        B_star = imf_complete(M_s, a, v_wind, t, Rs)
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
            P_rad = radio_power(P_in_mag, 0, nu, D)

            if IsBurst:
                I = I * (10 ** 1.53)
                P_rad *= 10**1.53

            intenss.append(P_rad)

    if flag:
        continue

    # print(f"{B_perp=}, {B_star=}, {n=}, {R_m=}. {Mdot=}, {P_in}, {P_rad}")
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

    observable_flag = False
    if 30 <= nu <= 180:
        for m in range(7):
            if Freq[m] <= nu <= Freq[m + 1]:
                if I >= (L_EU[m + 1] - L_EU[m]) / (Freq[m + 1] - Freq[m]) * (nu - Freq[m]) + L_EU[m]:
                    obs = str(EXO.name)
                    print(obs)
                    detectables.append(EXO)
                    observable_flag = True
                    break

    if not observable_flag:
        if 120 <= nu < 850 or 1050 < nu <= 1450:
            for m in range(4):
                if uGMRT["Frequencies"][m][0] <= nu <= uGMRT["Frequencies"][m][1] and I > uGMRT["RMS Noise"][m][0]:
                    obs = str(EXO.name)
                    print(obs)
                    detectables.append(EXO)
                    break

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
cond = np.where(semis < 2)

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

im = ax0.scatter(df1.x, df1.y, c=df1.s, s=df1.d, cmap="magma_r")
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
        ax0.annotate(txt, xy=(df1.x[i], df1.y[i]), xytext=(2, 2), textcoords="offset pixels", fontsize=7)

# lofar.plot(ax=ax0, x="Freq.", y="NL Core", style ="--", linewidth=0.2)
# lofar.plot(ax=ax0, x="Freq.", y="Full EU", style="g--", linewidth=0.5)
ax0.plot(Freq, L_EU, linestyle="dashed", color="purple", linewidth=0.5)

for i in range(4):
    x = uGMRT["Frequencies"][i]
    y = uGMRT["RMS Noise"][i]
    plt.plot(x, y, "b--", linewidth=0.5)
    if i == 0:
        ax0.fill_between(x, y, 10 ** 6, color="blue", alpha=0.1, label="uGMRT")
    else:
        ax0.fill_between(x, y, 10 ** 6, color="blue", alpha=0.1)

ax0.fill_between(Freq, L_EU, 10**6, color="purple", alpha=0.1, label="LOFAR")
plt.legend()
# ax0.fill_between(Freq, L_NL, 10**6, color="green", alpha=0.1)

fig0.colorbar(im, ax=ax0, label="Distance to Host Star ($\log_{10}{\mathrm{(AU)}}$)", aspect=25, extend="both")

ax0.axvline(x=10, color="black", linestyle="dashed")

ax0.set_xscale("log")
ax0.set_yscale("log")
# ax0.set_xlim(6, 30)
ax0.set_xlim(left=0.05)
ax0.set_ylim(bottom=10**(-10), top=10**0)

ax0.axvspan(0, 10, alpha=0.2, color="teal")
fig0.tight_layout()

if IsBurst:
    ax0.set_title("Frequency and Intensity of Burst CMI Emissions of the Exoplanet Sample")
else:
    ax0.set_title("Frequency and Intensity of Quiescent CMI Emissions of the Exoplanet Sample")
ax0.set_xlabel("Emission Frequency (MHz)")
ax0.set_ylabel("Radio Brightness (Jy)")

retro_noir(ax0)

TheBins1 = np.logspace(-12, 0, 13)

plt.rcParams['figure.figsize'] = [6, 4]
rc = {"font.size": 10}
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

TheBins2 = np.logspace(-3, 3, 13)

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

detectable_data = [[exo.name, exo.freq, exo.intensity] for exo in detectables]
with open("detectables.txt", "w") as fn:
    fn.write(tabulate(detectable_data))
    fn.close()

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
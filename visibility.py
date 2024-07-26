import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from adjustText import adjust_text
from radio_module import *
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun, get_moon
from astropy.time import Time
import astropy.units as u
import matplotlib.colors as mcolors
from datetime import datetime
# import smplotlib

plt.rcParams['figure.figsize'] = [7, 4]
# plt.rcParams['font.size'] = 10

df = pd.read_csv("obs_table.csv")

norm = mcolors.PowerNorm(gamma=0.3)  # Adjust gamma to control brightness

cmap = plt.cm.jet
cmap1 = plt.cm.turbo

which_data = np.load("observables.npz")
lofar = which_data["lofar"]
nenufar = which_data["nenufar"]
mwa = which_data["mwa"]
ugmrt = which_data["ugmrt"]

which = [lofar, nenufar, ugmrt, mwa]
tel_names = np.array(["LOFAR", "NenuFAR", "uGMRT", "MWA"])

LOFAR_lat = 54.9
NenuFAR_lat = 47.38
MWA_lat = -26.7
uGMRT_lat = 19.1

lats = np.array([LOFAR_lat, NenuFAR_lat, uGMRT_lat, MWA_lat])


def retro_noir(ax):
    ax.grid(alpha=0.1)
    ax.tick_params(direction="in", which="major", length=4, right=True, top=True, width=1.25)
    ax.tick_params(direction="in", which="minor", length=2, right=True, top=True, width=1.25)
    ax.spines['top'].set_linewidth(1.25)
    ax.spines['bottom'].set_linewidth(1.25)
    ax.spines['left'].set_linewidth(1.25)
    ax.spines['right'].set_linewidth(1.25)
    ax.minorticks_on()


def formatter(s):
    main, m, sec = s.split(":")
    val = abs(float(main)) + float(m)/60 + float(sec)/3600
    sign = s[0]
    if sign.isdigit() or sign == "+":
        return val
    elif sign == "-":
        return -val


def format_ra(ra_string):
    hours, minutes, seconds = ra_string.split()
    seconds = round(float(seconds))
    return f"{hours}:{minutes}:{seconds:02d}"


def two_format(x):
    return '{:.2f}'.format(x)


def three_format(x):
    return '{:.3f}'.format(x)


def one_format(x):
    return '{:.1f}'.format(x)


def zero_format(x):
    return '{:.0f}'.format(x)



all_df = df[["Name", "RA (J2000)", "DEC (J2000)"]].copy()

all_df["ra"] = all_df["RA (J2000)"].apply(formatter)
all_df["dec"] = all_df["DEC (J2000)"].apply(formatter)


def max_elev(ra, dec, obs_lat):
    return max(0, 90 - abs(dec - obs_lat))


def time_above_elevation(ra, dec, lat, min):
    dec, lat = dec/180*np.pi, lat/180*np.pi
    LST = np.linspace(0, 24, 1000)
    LHA = LST - ra
    LHA = LHA * 15 / 180 * np.pi
    h = np.arcsin(np.sin(lat) * np.sin(dec) + np.cos(lat) * np.cos(dec) * np.cos(LHA))
    cond = np.where(h > min/180 * np.pi)
    frac = len(h[cond]) / len(h)
    return frac*24


def possible_time_for_target(name, min=20, default=True, obs="LOFAR"):
    if default:
        tel = []
        latits = []
        times = []
        for ind, obs in enumerate(which):
            if name in obs:
                lat = lats[ind]
                latits.append(lat)
                tel.append(ind+1)
                ra, dec = all_df["ra"][df["Name"] == name], all_df["dec"][df["Name"] == name]
                t = time_above_elevation(ra, dec, lat, min)[0]
                times.append(t)
    else:
        tel = obs
        lat = lats[tel_names == tel]
    df_times = pd.DataFrame({"tel": tel,
                             "times": times})
    df_times = df_times.sort_values(by="times", ascending=False)
    tel_list = [x for x in df_times["tel"]]
    tel_list_str = ""
    for tel in tel_list:
        if tel_list_str != "" and tel_list.index(tel) != len(tel_list):
            tel_list_str += ", " + str(tel)
        else:
            tel_list_str += str(tel)

    latits = np.array(latits)
    lat = latits[np.where(np.array(times) == np.max(times))]
    ra, dec = all_df["ra"][df["Name"] == name], all_df["dec"][df["Name"] == name]
    h = time_above_elevation(ra, dec, lat, min)[0] // 1
    minutes = (time_above_elevation(ra, dec, lat, min)[0] - time_above_elevation(ra, dec, lat, min)[0] // 1)*60
    return tel_list_str, f"{h:.0f}h {minutes:.0f}m"


max_elev = np.vectorize(max_elev)
time_above_elevation = np.vectorize(time_above_elevation)
possible_time_for_target = np.vectorize(possible_time_for_target)
final = pd.read_csv("obs_table.csv")
final[["Telescope", "t20"]] = final["Name"].apply(lambda x: pd.Series(possible_time_for_target(x)))
final["a (AU)"] *= 1e2
final.rename(columns={"a (AU)": "a (10^-2 AU)"}, inplace=True)

# final["RA (J2000)"] = final["RA (J2000)"].apply(format_ra)
# final["DEC (J2000)"] = final["DEC (J2000)"].apply(format_ra)

final["Mass (MJ)"] = final["Mass (MJ)"].apply(three_format)
final["a (10^-2 AU)"] = final["a (10^-2 AU)"].apply(one_format)

final["Radius (RJ)"] = final["Radius (RJ)"].apply(two_format)
final["d (pc)"] = final["d (pc)"].apply(zero_format)
final["t (Gyr)"] = final["t (Gyr)"].apply(two_format)
final["vmax"] = final["vmax"].apply(zero_format)
final["Phi (mJy)"] = final["Phi (mJy)"].apply(two_format)


final.to_csv("final-table.csv", index=False)

# obs = [LOFAR_lat, NenuFAR_lat, MWA_lat, uGMRT_lat]
# obs_names = ["LOFAR", "NenuFAR", "MWA", "uGMRT"]

# obs = [NenuFAR_lat, LOFAR_lat, MWA_lat, uGMRT_lat]
# obs_names = ["NenuFAR", "LOFAR", "MWA", "uGMRT"]

obs = [LOFAR_lat, NenuFAR_lat, uGMRT_lat, MWA_lat]
obs_names = ["LOFAR", "NenuFAR", "uGMRT", "MWA"]

ras = np.linspace(0, 24, 200)
decs = np.linspace(-90, 90, 200)

ras, decs = np.meshgrid(ras, decs)


# plt.rcParams['figure.figsize'] = [9, 12]
#
# fig, axs = plt.subplots(len(obs), 2)
# ax_0, ax_1 = axs[:, 0], axs[:, 1]

results = []
for i in range(len(obs)):
    result = max_elev(ras, decs, obs[i])
    results.append(result)


def visibility_plot(results, save=False):

    plt.rcParams['font.size'] = 7  # 7 with NenuFAR
    fig, axs = plt.subplots(len(obs), 2, figsize=(7, 9))  # (7, 9) with NenuFAR, (9, 10) otherwise

    for i in range(len(results)):

        result = results[i]

        ax_0, ax_1 = axs[:, 0], axs[:, 1]

        ax = ax_0[i]
        ax1 = ax_1[i]

        # cp = ax.contourf(X, Y, L_r)
        cp = ax.imshow(result, cmap=cmap, origin="lower",
                        extent=[ras.min(), ras.max(), decs.min(), decs.max()], aspect="auto")
        fig.colorbar(cp, label=r"Max. Elevation Angle ($\degree$)", ax=ax, extend="both")  # Add a colorbar to a plot
        ax.set_title(obs_names[i])
        ax.set_xlabel(r'R.A. (h)')
        ax.set_ylabel(r'DEC ($\degree$)')

        df1 = all_df[all_df["Name"].isin(which[i])].reset_index(drop=True)
        for index, row in df1.iterrows():
            # Extract the values from the current row
            name = row['Name']
            ra = row['ra']
            dec = row['dec']

            # Add a circular point at coordinates (ra, dec) with white face color and black edge color
            ax.plot(ra, dec, 'o', color='white', markersize=7, markeredgewidth=1, markeredgecolor='black')
            ax1.plot(ra, dec, 'o', color='white', markersize=7, markeredgewidth=1, markeredgecolor='black')

            # Annotate the point with the name above the circular point in white color
            # ax.annotate(name, (ra, dec), xytext=(0, 5), textcoords='offset points', ha='center', color='white')

        texts = [ax.text(df1.ra[k], df1.dec[k], df1.Name[k], ha='center', va='center', color="white") for k in range(len(df1["Name"]))]
        texts1 = [ax1.text(df1.ra[k], df1.dec[k], df1.Name[k], ha='center', va='center', color="white") for k in range(len(df1["Name"]))]

        res1 = time_above_elevation(ras, decs, obs[i], 20)
        cp = ax1.imshow(res1, cmap=cmap1, origin="lower",
                       extent=[ras.min(), ras.max(), decs.min(), decs.max()], aspect="auto")
        fig.colorbar(cp, label=r"Time above $20\degree$ (Hours/day)", ax=ax1, extend="both")  # Add a colorbar to a plot
        ax1.set_title(obs_names[i])
        ax1.set_xlabel(r'R.A. (h)')
        ax1.set_ylabel(r'DEC ($\degree$)')

        retro_noir(ax)
        retro_noir(ax1)

        adjust_text(texts, ax=ax)
        adjust_text(texts1, ax=ax1)

    plt.subplots_adjust(hspace=0.6)

    plt.show()

    if save:
        plt.savefig("visibility.pdf")

# time_above_elevation(3.1, 17.2, 39.5, 0)


visibility_plot(results)

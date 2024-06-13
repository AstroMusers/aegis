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


df1 = df[["Name", "RA (J2000)", "DEC (J2000)"]].copy()

df1["ra"] = df1["RA (J2000)"].apply(formatter)
df1["dec"] = df1["DEC (J2000)"].apply(formatter)


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



LOFAR_lat = 54.9
MWA_lat = -26.7
uGMRT_lat = 19.1

obs = [LOFAR_lat, MWA_lat, uGMRT_lat]
obs_names = ["LOFAR", "MWA", "uGMRT"]

ras = np.linspace(0, 24, 200)
decs = np.linspace(-90, 90, 200)

ras, decs = np.meshgrid(ras, decs)

elev = np.vectorize(max_elev)
time_above_elevation = np.vectorize(time_above_elevation)

plt.rcParams['figure.figsize'] = [10, 15]

fig, axs = plt.subplots(len(obs), 2)
ax_0, ax_1 = axs[:, 0], axs[:, 1]

for i in range(len(obs)):

    results = elev(ras, decs, obs[i])

    ax = ax_0[i]
    ax1 = ax_1[i]

    # cp = ax.contourf(X, Y, L_r)
    cp = ax.imshow(results, cmap=cmap, origin="lower",
                    extent=[ras.min(), ras.max(), decs.min(), decs.max()], aspect="auto")
    fig.colorbar(cp, label=r"Max. Elevation Angle ($\degree$)", ax=ax, extend="both")  # Add a colorbar to a plot
    ax.set_title(obs_names[i])
    ax.set_xlabel(r'R.A. (h)')
    ax.set_ylabel(r'DEC ($\degree$)')

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

    texts = [ax.text(df1.ra[k], df1.dec[k], df1.Name[k], ha='center', va='center', fontsize=9, color="white") for k in range(len(df1["Name"]))]
    texts1 = [ax1.text(df1.ra[k], df1.dec[k], df1.Name[k], ha='center', va='center', fontsize=9, color="white") for k in range(len(df1["Name"]))]

    res1 = time_above_elevation(ras, decs, obs[i], 20)
    cp = ax1.imshow(res1, cmap=cmap1, origin="lower",
                   extent=[ras.min(), ras.max(), decs.min(), decs.max()], aspect="auto")
    fig.colorbar(cp, label=r"Time above $20\degree$ (Hours/day)", ax=ax1, extend="both")  # Add a colorbar to a plot
    ax1.set_title(obs_names[i])
    ax1.set_xlabel(r'R.A. (h)')
    ax1.set_ylabel(r'DEC ($\degree$)')

    plt.subplots_adjust(hspace=0.6)

    retro_noir(ax)
    retro_noir(ax1)

    adjust_text(texts, ax=ax)
    adjust_text(texts1, ax=ax1)


# time_above_elevation(3.1, 17.2, 39.5, 0)

plt.show()

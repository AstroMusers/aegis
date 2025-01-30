import numpy as np
import pandas as pd
import math
from adjustText import adjust_text
from radio_module import *
from rotation_script import *
from copy import deepcopy
import smplotlib
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker
import json

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

IsBurst = 1

df = pd.read_csv("df1.csv")
df.fillna('', inplace=True)
all = np.load("all.npz", allow_pickle=True)

y_both_maxerr, y_both_minerr, y_kin_minerr, y_kin_maxerr, y_mag_minerr, y_mag_maxerr = all["y_both_maxerr"], all["y_both_minerr"], all["y_kin_minerr"], all["y_kin_maxerr"], all["y_mag_minerr"], all["y_mag_maxerr"]
x_minerr, x_maxerr = all["x_minerr"], all["x_maxerr"]
y_mag_err = [y_mag_minerr, y_mag_maxerr]
y_kin_err = [y_kin_minerr, y_kin_maxerr]
y_both_err = [y_both_minerr, y_both_maxerr]
y_err = [y_mag_err, y_kin_err, y_both_err]
x_err = [x_minerr, x_maxerr]

detectables_both = list(all["detectables_both"])
average_errors = list(all["average_errors"])
intensities, magnetic_fields = list(all["intensities"]), list(all["magnetic_fields"])

real_outliers = set(all["real_outliers"])

det_data = [[exo.name, exo.freq, exo.intensity_both] for exo in detectables_both]
df_det = pd.DataFrame(det_data[0:], columns=["Name", "Freq", "Flux"])

def minor_tick_format(x, pos):
    if x in [i * 10**j for j in range(-1, 3) for i in range(2, 10, 2)]:  # Customize range as needed
        return f"{x:g}"  # Format in plain numbers
    return ""

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

    plt.show()


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
                                           c=non_outliers.s, cmap='magma_r', marker='o',
                                           norm=Normalize(vmin=min_color, vmax=max_color))

        # Plot outliers with upside-down triangles, using size and color mappings
        scatter_outliers = ax0.scatter(outliers.xerr_nonzero, outliers.yerr_nonzero, s=outliers.d,
                                       c=outliers.s, cmap='magma_r', marker='v',
                                       norm=Normalize(vmin=min_color, vmax=max_color))

    else:
        scatter_non_outliers = ax0.scatter(non_outliers.xerr_nonzero, non_outliers.yerr_nonzero,
                                           s=size[~df['Names'].isin(real_outliers)],
                                           c=non_outliers.s, cmap='magma_r', marker='o',
                                           norm=Normalize(vmin=min_color, vmax=max_color))

        # Plot outliers with upside-down triangles, using size and color mappings
        scatter_outliers = ax0.scatter(outliers.xerr_nonzero, outliers.yerr_nonzero,
                                       s=size[df['Names'].isin(real_outliers)],
                                       c=outliers.s, cmap='magma_r', marker='v',
                                       norm=Normalize(vmin=min_color, vmax=max_color))

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

    # x_nenu = np.linspace(NenuFreq[0], NenuFreq[1], 100)
    # y_nenu = np.linspace(NenuNoise[0], NenuNoise[1], 100)
    x_nenu = NenuFreq
    y_nenu = NenuNoise
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
    fig0.colorbar(sm, ax=ax0, label=r"$\log_{10}$" + "(Semi‐major Axis [AU])", aspect=25, extend="both")

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
            texts = [plt.text(df.x[i], df.y[i], lab_strict[i], ha='center', va='center', fontsize=8) for i in
                     range(len(lab_strict)) if lab_strict[i] != ""]
        elif fix_lim:
            texts = [plt.text(df.x[i], df.y[i], lab[i], ha='center', va='center', fontsize=8) for i in range(len(lab))
                     if lab[i] != "" and is_within_limits(df.x[i], df.y[i], xlim, ylim)]
        else:
            texts = [ax0.text(df.x[i]*1.13, df.y[i]*1.10, lab[i], ha='center', va='center', fontsize=8) for i in range(len(lab))
                     if lab[i] != ""]
            ax0.xaxis.set_minor_formatter(ticker.FuncFormatter(minor_tick_format))
            ax0.xaxis.set_major_formatter(ticker.LogFormatter())
        fig0.legend(fontsize=13, bbox_to_anchor=(0.1, 0.15), loc="lower left", frameon=True)
        line1 = Line2D([0], [0], marker="v", linestyle="None", markerfacecolor="orange", markeredgecolor="black")
        line2 = Line2D([0], [0], marker="o", linestyle="None", markerfacecolor="orange", markeredgecolor="black")
        # fig0.legend((line1, line2), ("Outliers", "Insiders"), frameon=True, shadow=True, bbox_to_anchor=(0.8, 0.95), fontsize=12)
        fig0.legend((line1,), ("Outliers",), frameon=True, shadow=True, bbox_to_anchor=(0.8, 0.95), fontsize=12)

        # adjust_text(texts, arrowprops=dict(arrowstyle="-", color="k", lw=0.5),
        #             force_points=(3, 3), force_text=(2, 2), force_objects=(1.5, 1.5),
        #             expand_points=(1.15, 1.15), expand_objects=(1.5, 1.5), expand_align=(1.2, 1.2), precision=20)

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
        ax0.annotate("", xy=(center_x, center_y * avg_err[2]), xytext=(center_x, center_y * avg_err[3]),
                     arrowprops=arrowprops)
        ax0.annotate("", xy=(center_x * avg_err[0], center_y), xytext=(center_x * avg_err[1], center_y),
                     arrowprops=arrowprops)

        ax0.legend(loc="center right", bbox_to_anchor=(1, 0.3), fontsize=11, frameon=True, shadow=True, ncol=1)

        xmin, xmax = ax0.get_xlim()

        # ax0.annotate("Observable From Ground", xy=(np.sqrt(xmax*10), max(df.y)*1.5),
        #     xytext=(np.sqrt(xmax*10), max(df.y)*2), ha='center',  # Horizontal alignment of the text
        #     fontsize=12,  # Size of the text
        #              )

        # ax0.annotate("", xy=(10, max(df.y)*1.5), xytext=(xmax, max(df.y)*1.5), arrowprops=dict(arrowstyle='<-', color='black', linestyle="dotted"))

        ax0.text(np.sqrt(xmax * 10), min(df.y) / 4, "Observable From Ground", color="green",
                 ha='center', va="center",  # Horizontal alignment of the text
                 fontsize=12,  # Size of the text
                 bbox=dict(facecolor='none', edgecolor='green', boxstyle='square')
                 )

        ax0.text(np.sqrt(xmin * 10), min(df.y) / 4, "Cannot Penetrate the Ionosphere", color="red",
                 ha='center', va="center",  # Horizontal alignment of the text
                 fontsize=12,  # Size of the text
                 bbox=dict(facecolor='none', edgecolor='red', boxstyle='sawtooth')
                 )

    ax0.set_xlabel("Maximum Emission Frequency [MHz]")
    ax0.set_ylabel("Radio Flux Density [Jy]")
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

    plt.show()

scatter_plot(df, "both", y_err, x_err, df_det, average_errors, save=True)
scatter_plot(df, "both", y_err, x_err, df_det, average_errors, save=True, zoom=True, others=10)
outcome_dist_hists(intensities, "both", magnetic_fields)

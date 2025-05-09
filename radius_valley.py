import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import smplotlib

filename = "NASA1904.csv"
df = pd.read_csv(filename, comment="#")
df2 = pd.read_csv("Output Tables/all.csv")
reachers = df[df["pl_name"].isin(df2["Name"])]
enough_data = reachers[["pl_name", "pl_bmassj", "pl_radj", "pl_orbsmax", "sy_dist", "st_age", "pl_orbper"]]
enough_data.rename(columns={"pl_name": "Name"}, inplace=True)
merged = pd.merge(enough_data, df2, on="Name")

merged_close = merged[(merged["pl_orbsmax"] < 0.1) & (merged["pl_radj"] * 11.2089 < 7)]
# merged_close = merged[(merged["pl_orbsmax"] < 0.1) & (merged["pl_orbper"] < 100)]

radii = np.array(merged_close["pl_radj"]) * 11.2089
bFields = np.array(merged_close["Max. Frequency [MHz]"] * 14 / 40)
masses = np.array(merged_close["pl_bmassj"])
semis = np.array(merged_close["pl_orbsmax"])
pers = np.array(merged_close["pl_orbper"])
ages = np.array(merged_close["st_age"])

irradiance = np.array(7.02e-2 * (merged_close["st_age"]/4.6) ** (-1.74) / merged_close["pl_orbsmax"]**2)

x = bFields
y = radii

group_masker = irradiance
masks = [np.where(group_masker < 10), np.where((group_masker > 10) & (group_masker < 100)), np.where(group_masker > 100)]

color1 = "xkcd:deep red"
color2 = "xkcd:indigo"

# Create the figure and the subplots
# plt.rcParams['font.size'] = 21
fig, axs = plt.subplots(len(masks), 2, figsize=(5.25, 7.5), gridspec_kw={'width_ratios': [1, 1]}, sharey=True, sharex="col")
groups = [r"$R/R_\oplus$(Low $F_X$)", r"$R/R_\oplus$(Moderate $F_x$)", r"$R/R_\oplus$(High $F_X$)"]

for i, mask in enumerate(masks):
    masker = pers[mask]
    median_masker = np.median(masker)
    colors_masked = np.where(masker > median_masker, color1, color2)

    ax1, ax2 = axs[i][0], axs[i][1]
    x_masked = x[mask]
    y_masked = y[mask]

    # Scatter plot on the left
    for color in [color1, color2]:
        idx = np.where(colors_masked == color)
        ax1.scatter(x_masked[idx], y_masked[idx], marker="o", linewidth=0.7, facecolors="none", edgecolors=color)

    ax1.set_xscale('log')
    ax1.set_ylabel(groups[i])
    ax1.tick_params(which="major", length=12)
    ax1.tick_params(which="minor", length=6)
    ax1.grid("on", alpha=0.2)

    # Histogram on the right
    bins = np.linspace(0, np.max(y_masked) // 1 + 1, 25)
    # bins = np.linspace(0, np.max(y_masked) // 1 + 1, 26)
    bins = np.linspace(0, 7, 25)
    # bins = np.arange(0, 7.25, 0.25)
    label_long = r"Longer $T_\mathrm{orb}$"
    label_short = r"Shorter $T_\mathrm{orb}$"
    ax2.hist(y_masked[colors_masked == color1], bins=bins, orientation='horizontal', color=color1, alpha=0.8, edgecolor=color1, linewidth=1, histtype="step", hatch="//////", label=label_long)
    ax2.hist(y_masked[colors_masked == color2], bins=bins, orientation='horizontal', color=color2, alpha=0.8, edgecolor=color2, linewidth=1, histtype="step", hatch="\\\\\\\\\\\\", label=label_short)
    ax2.hist(y_masked, bins=bins, orientation='horizontal', color='k', edgecolor='black', linewidth=2, histtype="step")
    ax2.tick_params(which="major", length=15)
    ax2.tick_params(which="minor", length=8)

    if i == len(masks)-1:
        ax1.set_xlabel(r'$B_\mathrm{surf}$ [G]')
        ax2.set_xlabel('N(R)')
    ax2.yaxis.tick_right()

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.05, wspace=0.05)
    if i == 1:
        ax2.legend(loc=1, fontsize=12)
plt.savefig("valleys.pdf")
plt.show()


# Begin single valley plot that does not go into the paper

colors = np.where(pers > np.median(pers), color1, color2)
plt.rcParams['font.size'] = 19

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7, 5), gridspec_kw={'width_ratios': [1, 1]}, sharey=True)
# Scatter plot on the left
for color in [color1, color2]:
    idx = np.where(colors == color)
    ax1.scatter(x[idx], radii[idx], marker="o", linewidth=0.7, facecolors="none", edgecolors=color)

ax1.set_xscale('log')
ax1.set_xlabel(r'$B_\mathrm{surf}$ [G]')
ax1.set_ylabel(r'Radius $[R_\oplus]$')
ax1.tick_params(which="major", length=12)
ax1.tick_params(which="minor", length=6)
ax1.grid("on", alpha=0.2)

# Histogram on the right
bins = np.linspace(0, 7, 25)
ax2.hist(radii[colors == color1], bins=bins, orientation='horizontal', color=color1, alpha=0.8, edgecolor=color1, linewidth=1, histtype="step", hatch="//////", label=r"Longer $T_\mathrm{orb}$")
ax2.hist(radii[colors == color2], bins=bins, orientation='horizontal', color=color2, alpha=0.8, edgecolor=color2, linewidth=1, histtype="step", hatch="\\\\\\\\\\\\", label=r"Shorter $T_\mathrm{orb}$")
ax2.hist(radii, bins=bins, orientation='horizontal', color='k', edgecolor='black', linewidth=2, histtype="step")
ax2.tick_params(which="major", length=15)
ax2.tick_params(which="minor", length=8)

ax2.set_xlabel('N(R)')
ax2.yaxis.tick_right()

plt.tight_layout()
plt.legend(fontsize=12)
plt.savefig("single_valley.pdf")
plt.show()


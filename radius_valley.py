import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import smplotlib

filename = "NASA1407.csv"
df = pd.read_csv(filename, comment="#")
df2 = pd.read_csv("Output Tables/all.csv")
reachers = df[df["pl_name"].isin(df2["Name"])]
enough_data = reachers[["pl_name", "pl_bmassj", "pl_radj", "pl_orbsmax", "sy_dist", "st_age", "pl_orbper"]]
enough_data.rename(columns={"pl_name": "Name"}, inplace=True)
merged = pd.merge(enough_data, df2, on="Name")

merged_close = merged[(merged["pl_orbsmax"] < 100) & (merged["pl_radj"] * 11.2089 < 7)]
# merged_close = merged[(merged["pl_orbsmax"] < 0.1) & (merged["pl_orbper"] < 100)]

radii = np.array(merged_close["pl_radj"]) * 11.2089
bFields = np.array(merged_close["Max. Frequency [MHz]"] * 14 / 40)
masses = np.array(merged_close["pl_bmassj"])
semis = np.array(merged_close["pl_orbsmax"])
pers = np.array(merged_close["pl_orbper"])

x = bFields
y = radii

# Determine the median of the flux values
median_x = np.median(x)

color1 = "xkcd:deep red"
color2 = "xkcd:indigo"

# Assign colors based on whether the flux is above or below the median
colors = np.where(pers > np.median(pers), color1, color2)
# colors = np.where(semis < 0.1, color1, color2)
# median_constrained = np.median(x)
# colors = np.where(x > 80, color1, color2)


# Create the figure and the subplots
plt.rcParams['font.size'] = 21
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
bins = np.linspace(0, np.max(radii) // 1 + 1, 25)
# bins = np.linspace(0, 7.5, 35)
ax2.hist(radii[colors == color1], bins=bins, orientation='horizontal', color=color1, alpha=0.8, edgecolor=color1, linewidth=1, histtype="step", hatch="//////", label=r"Longer $T_\mathrm{orb}$")
ax2.hist(radii[colors == color2], bins=bins, orientation='horizontal', color=color2, alpha=0.8, edgecolor=color2, linewidth=1, histtype="step", hatch="\\\\\\\\\\\\", label=r"Shorter $T_\mathrm{orb}$")
ax2.hist(radii, bins=bins, orientation='horizontal', color='k', edgecolor='black', linewidth=2, histtype="step")
ax2.tick_params(which="major", length=15)
ax2.tick_params(which="minor", length=8)

ax2.set_xlabel('N(R)')
# ax2.set_ylabel('Radius $[R_\oplus]$')
# ax2.set_ylim(0, 7.5)
# ax2.yaxis.set_label_position("right")
ax2.yaxis.tick_right()

plt.tight_layout()
plt.legend(fontsize=17)
plt.show()

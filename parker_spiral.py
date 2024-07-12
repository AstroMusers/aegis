import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib.lines import Line2D
import smplotlib
import pandas as pd
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

data = np.load("taub_wind.npz")

ranges = data["ranges"]
speeds = data["speeds"]

colormap = "viridis"

M = 1.4  # Solar Mass
R = 1.44  # Solar Radius
G = 8.87129 * 10**2  # Gravitational constant in AU (M_sun)-1 (km/s)2
t = 2  # Gyr
coeff = 1.89 / (4.6**(-0.655))
B0 = coeff * t**(-0.655)  # Gauss
P_rot = 4  # Days

a = 0.049  # AU

conversions_constant = 2 * np.pi * 1731  #

ranges = ranges.tolist()
speeds = speeds.tolist()

ranges = np.array(ranges)
speeds = np.array(speeds)

kepler = np.sqrt(G*M/ranges)
v_eff = np.sqrt(speeds**2 + kepler**2),

theta = np.linspace(0, 2 * np.pi, 101)
r, theta = np.meshgrid(ranges, theta)

B_r = B0 * ((R/215) / r)**2

B_phi = B_r * conversions_constant * r / P_rot / speeds

B = np.sqrt(B_r**2 + B_phi**2)

alpha = np.arctan(B_phi / B_r)
beta = np.arctan(kepler / speeds)

B_perp = B * abs(np.sin(alpha - beta)) * 1e5  # nT

B_perp = np.log10(B_perp)
B = np.log10(B) + 5  # nT

x = r * np.cos(theta)
y = r * np.sin(theta)

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, constrained_layout=True)
ax.grid(False, which='minor')

r = np.log10(r)

masked_values = np.ma.masked_where(r < np.log10(R/215), B_perp)

# Plot the data
c = ax.pcolormesh(theta, r, masked_values, shading='gouraud', cmap=colormap)
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

ax.grid("on", color="w", alpha=0.25,)

# Create a new axis for the colorbar
# cbar_ax = fig.add_axes([0.1, 0, 0.8, 0.03])  # [left, bottom, width, height]

# Add the colorbar to the new axis
# fig.colorbar(c, cax=cbar_ax, orientation="horizontal", extend="both")
fig.colorbar(c, orientation="vertical", pad=0.09, extend="both", fraction=0.03, aspect=30, label="log$_{10}$($B_\perp$ [nT])")

# Set labels
ax.set_xlabel('$\phi$', loc="right", labelpad=-180, color="white")
# ax.set_ylabel('log(Radius (AU))', labelpad=50)

ax.set_rlabel_position(-90)  # Move the radial labels away from the plotted data

ax.tick_params(axis='y', colors='pink')  # Radial ticks
ax.set_rticks([-1, -0.5, 0])

# Optional: Label the radial axis with context-specific information
ax.text(- 95 * np.pi / 180, (np.min(r)+np.max(r))/2 + 0.8, 'log$_{10}$($r$ [AU])', ha='center', va='center', fontsize=13, rotation=90, color="pink")
# ax.text(np.pi, (np.min(r)+np.max(r))/2 - 0.8, 'log($r$ [AU])', ha='center', va='center', fontsize=13)

circle_radius = -r[0][0] + np.log10(a)
circle = Circle((0, 0), circle_radius, transform=ax.transData._b, color='w', fill=False, linestyle='-.', linewidth=2)
ax.add_patch(circle)

circle2_radius = -r[0][0] + np.log10(R/215)
circle2 = Circle((0, 0), circle2_radius, transform=ax.transData._b, color='red', fill=False, linestyle='-', linewidth=2)
ax.add_patch(circle2)

dot_radius = np.log10(a)
dot_angle = 120 * np.pi / 180  # Example angle for the dot position
dot_x = dot_radius * np.cos(dot_angle)
dot_y = dot_radius * np.sin(dot_angle)
ax.scatter(dot_angle, dot_radius, color="red", edgecolor="k")  # 'ro' specifies a red dot
# ax.plot(0, np.min(r), "*", color="orange", markersize=15)
ax.text(dot_angle, dot_radius+0.2, "$\\tau$ Boo b", ha="center", va="center", color="yellow", rotation=dot_angle * 180 / np.pi - 90)

line1 = Line2D([], [], color="w", linestyle='-.', linewidth=2)
line2 = Line2D([], [], color="red", linestyle='-', linewidth=2)


def add_curved_arrow(ax, radius, start_angle, end_angle):
    arrow = FancyArrowPatch((start_angle, radius), (end_angle, radius),
                            transform=ax.transData, color='white', connectionstyle=f"arc3,rad=0.12",
                            arrowstyle='->', mutation_scale=15)
    ax.add_patch(arrow)


# Example arrows
add_curved_arrow(ax, 0.1, 0, 15 * np.pi / 180)
# add_curved_arrow(ax, circle_radius, 3 * np.pi / 2, 7 * np.pi / 4)

from scipy.interpolate import griddata

# Flatten the meshgrid data
points = np.array([r.flatten(), theta.flatten()]).T
values_flat = masked_values.flatten()

# Point where we want to find the color
point_of_interest = np.array([[dot_radius, dot_angle]])

# Interpolate the value at this point
interpolated_value = griddata(points, values_flat, point_of_interest, method='linear')

# Get the color corresponding to this interpolated value
norm = Normalize(vmin=masked_values.min(), vmax=masked_values.max())
cmap = plt.get_cmap(colormap)
mappable = ScalarMappable(norm=norm, cmap=cmap)
dot_color = mappable.to_rgba(interpolated_value[0])


# fig.legend((line1, line2), ("Orbit", "Stellar Surface"), frameon=True, shadow=True, facecolor=dot_color, labelcolor="white", ncol=2)
fig.legend((line1,), ("Orbit",), frameon=True, shadow=True, facecolor=dot_color, labelcolor="w", ncol=1)
fig.legend((line2,), ("Stellar Surface",), frameon=True, shadow=True, labelcolor="k", ncol=1, loc=2)

# fig.tight_layout()
fig.savefig("parker_spiral.pdf", dpi=300)

plt.show()


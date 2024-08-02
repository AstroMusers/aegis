import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Arc
import matplotlib as mpl
import smplotlib


# Note! By convention, polar angle is measured from the equator instead of the poles. Caution is required when using mainstream physics conventions

# Function to distort the dipole field line under the influence of stellar wind. Just an approximation.

factor = 0.7


def distorted_field_line(theta, r_0, r_s):
    # Compression on the dayside (theta in [0, pi]) and elongation on the nightside (theta in [pi, 2*pi])
    f_theta = 1 + factor * np.cos(theta)
    return r_0 * np.cos(theta)**2 * f_theta


# Example R_mp
r_s = 3.84

# Define theta array
theta = np.linspace(0, 2 * np.pi, 2000)

# Generate field lines for different values of r_0
r_0_values = np.linspace(r_s / 100, r_s, 10)

# plt.rcParams['font.size'] = 21
plt.rcParams["font.size"] = 18
fig, ax = plt.subplots(figsize=(7, 5))

planet_radius = 1

# Plot distorted field lines
xs, ys = [], []
for i, r_0 in enumerate(r_0_values):
    r = distorted_field_line(theta, r_0, r_s) / (1-factor)
    mask = np.where(r > planet_radius*1.15)
    r_masked = r[mask]
    theta_masked = theta[mask]
    x = r_masked * np.cos(theta_masked)
    y = r_masked * np.sin(theta_masked)
    xs.append(x)
    ys.append(y)

    if i == len(r_0_values) - 1:
        ax.scatter(x, y, color="k", marker="_", s=2)
    else:
        ax.scatter(x, y, color="k", marker="_", s=2)

# test = np.min(np.array(xs))

planet = Circle((0, 0), planet_radius, facecolor='xkcd:sea', edgecolor="k", alpha=0.9)
ax.add_patch(planet)

r_s = -r_s
# Plot the standoff distance
# plt.plot([-r_s * 1e-7, -r_s* 1e-7], [-2, 2], 'k--', label=f'Standoff distance = {r_s:.2f} R_s')

lam = np.arccos(np.sqrt(1/-r_s))
x_em, y_em = planet_radius * np.cos(lam), planet_radius * np.sin(lam)


# R_mp
start = (0, 0)
end = (r_s, 0)
ax.annotate("", xy=end, xytext=start, arrowprops=dict(arrowstyle='->', lw=1.5, color="red"))
ax.annotate(r"$R_\mathrm{mp}$", xy=(0, 0.1), xytext=(r_s/2-0.05, 0.1), color="r", fontsize=14)

# Arrow to lambda
ax.annotate("", xy=(-x_em, y_em), xytext=(0, 0), arrowprops=dict(arrowstyle='->', lw=1.5, color="red"))
ax.annotate(r"$R$", xy=(0, 0.1), xytext=(-x_em/2, y_em/2), color="r", fontsize=14)

# Lambda
ax.annotate("",  xy=(x_em, y_em), xytext=(-x_em, y_em), arrowprops=dict(arrowstyle='-', lw=1.5, color="yellow"))
ax.annotate("",  xy=(x_em, -y_em), xytext=(-x_em, -y_em), arrowprops=dict(arrowstyle='-', lw=1.5, color="yellow"))
center_x, center_y = 0, 0
width, height = 1, 1
angle = 0  # Rotation angle of the arc
end_angle = 180  # Starting angle of the arc
start_angle = 180 - lam*180/np.pi  # Ending angle of the arc

# Create the arc
arc = Arc((center_x, center_y), width, height, angle=angle, theta1=start_angle, theta2=end_angle, edgecolor='yellow', lw=0.5)
ax.add_patch(arc)
ax.annotate(r"$\lambda_\mathrm{CMI}$", xy=(0, 0), xytext=(-0.7, 0.1), color="yellow", fontsize=14)

# Wind
x_wind = np.linspace(-6, -4.5, 100)
y_for_wind = np.linspace(-4, 4, 10)

for i, y in enumerate(y_for_wind):
    y_wind = np.linspace(y, y, 100)

    with mpl.rc_context({'path.sketch': (3, 15, 1)}):
        ax.plot(x_wind, y_wind, color="xkcd:orange", lw=2)
    ax.annotate("", xy=(x_wind[-1]*0.88, y_wind[-1]), xytext=(x_wind[-1], y_wind[-1]), arrowprops=dict(arrowstyle='->', lw=2, color="xkcd:orange"))

ax.annotate(r"$v_\mathrm{sw}$", xy=(0.15, 0.9), xytext=((x_wind[0] + x_wind[-1])/1.87, y+0.2), color="xkcd:orange", fontsize=16)


# Calculate beaming angle associated with a band of emission on the planet where the total solid angle is 1.6 sr.
solid_angle = 1.6 / 2
beam_angle = np.arcsin(solid_angle / (4*np.pi*np.cos(lam)))

rays = np.linspace(1, 6, 1000)
angles = np.linspace(lam-beam_angle, lam+beam_angle, 500)
combinations = [(-1, -1), (-1, 1), (1, -1), (1, 1)]
for angle in angles:
    for couple in combinations:
        c1, c2 = couple[0], couple[1]
        xs = rays * np.cos(angle) * c1
        ys = rays * np.sin(angle) * c2
        ax.plot(xs, ys, color="purple", alpha=0.02)


ax.set(xlabel="X [R]", ylabel="Y [R]")
text_str = r"$R_\mathrm{mp} = $" + f"{-r_s:.2f}" + r"$\,R$"
props = dict(boxstyle='square', facecolor='none', alpha=1)

# Using coordinates (1, 1) for upper right corner with offset
ax.text(0.5, 0.95, text_str, transform=ax.transAxes, fontsize=18,
        verticalalignment='top', horizontalalignment='center', bbox=props)

ax.set(xlim=[r_s*1.75, -r_s*1.75], ylim=[r_s*1.25, -r_s*1.25])


plt.show()

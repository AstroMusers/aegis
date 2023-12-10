import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

rc = {"font.family": "times new roman",
      "font.size": 11,
      "mathtext.fontset": "stix"}
plt.rcParams.update(rc)

mpl.rcParams["figure.autolayout"] = True
plt.rcParams['figure.figsize'] = [4, 9]

mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=["navy", "darkred", "b"])

R_j = 69911  # Jupiter radius in kilometers
mu_0 = 4 * np.pi * 10 ** (-7)  # Vacuum Permeability in SI units

plt.figure(1)

v_sw = np.linspace(1, 1000, 1000)  # km/s
R_mp = 56  # Jupiter Radius
B_sw = 10 ** (-8)

L_r = 3 * 10 ** (-3) * np.pi * (R_mp * R_j * 10 ** 3) ** 2 * (v_sw * 10 ** 3) * B_sw ** 2 / (2 * mu_0)

plt.subplot(3, 1, 1)
plt.plot(v_sw, L_r)
plt.xlabel(r"$v_{sw}$ (km/s)")
plt.ylabel(r"$L_r$ (Watts)")
plt.title(r"Scaling of $L_r$ by $v_{sw}$")

v_sw = 400  # km/s
R_mp = np.linspace(0.1, 100, 1000)  # Jupiter radius
B_sw = 10 ** (-8)

L_r = 3 * 10 ** (-3) * np.pi * (R_mp * R_j * 10 ** 3) ** 2 * (v_sw * 10 ** 3) * B_sw ** 2 / (2 * mu_0)

plt.subplot(3, 1, 2)
plt.plot(R_mp, L_r)
plt.xlabel(r"$R_{mp}$ $(R_J)$")
plt.ylabel(r"$L_r$ (Watts)")
plt.title(r"Scaling of $L_r$ by $R_{mp}$")

v_sw = 400  # km/s
R_mp = 56  # Jupiter Radius
B_sw = np.linspace(1, 100, 1000)  # nT

L_r = 3 * 10 ** (-3) * np.pi * (R_mp * R_j * 10 ** 3) ** 2 * (v_sw * 10 ** 3) * (B_sw * 10 ** (-9)) ** 2 / (2 * mu_0)

plt.subplot(3, 1, 3)
plt.plot(B_sw, L_r)
plt.xlabel(r"$B_{sw}$ (nT)")
plt.ylabel(r"$L_r$ (Watts)")
plt.title(r"Scaling of $L_r$ by $B_{sw}$")

plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))


mpl.rcParams["figure.autolayout"] = True
plt.rcParams['figure.figsize'] = [5, 3]

rc = {"font.family": "times new roman",
      "font.size": 11,
      "mathtext.fontset": "stix"}
plt.rcParams.update(rc)

mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=["navy", "darkred", "b"])

k = 7.39 * 10 ** 22

v_sw = 400  # km/s
B_sw = 10 ** (-8)  # teslas
B_s = np.linspace(1, 50, 1000)
a = np.linspace(0.1, 10, 1000)

B_s, a = np.meshgrid(B_s, a)

L_r = k * B_s ** (2 / 3) * a ** (2 / 3) * v_sw * B_sw ** 2

fig, ax = plt.subplots(1, 1)
# cp = ax.contourf(X, Y, L_r)
cp = plt.imshow(L_r, cmap=plt.cm.viridis, origin="lower",
                extent=[B_s.min(), B_s.max(), a.min(), a.max()], aspect="auto")
fig.colorbar(cp, label=r"$L_r$ (Watts)")  # Add a colorbar to a plot
ax.set_title('Scaling of $L_r$ by $B_s$ and $a$')
ax.set_xlabel(r'$B_s$ (G)')
ax.set_ylabel(r'$a$ (AU)')
plt.show()

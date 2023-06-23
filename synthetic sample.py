import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import pandas as pd

mpl.rcParams["figure.autolayout"] = True
plt.rcParams['figure.figsize'] = [8, 4]

rc = {"font.family": "times new roman",
      "font.size": 11,
      "mathtext.fontset": "stix"}
plt.rcParams.update(rc)

mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=["firebrick", "navy", "b"])

rng = np.random.default_rng()


class Exoplanet:
    def __init__(self, semi_major_axis, magnetic_field, star_mass, star_mass_loss, distance):
        self.semi_major_axis = semi_major_axis
        self.magnetic_field = magnetic_field
        self.star_mass = star_mass
        self.star_mass_loss = star_mass_loss
        self.distance = distance

    def __repr__(self):
        return f"Exoplanet with: a={self.semi_major_axis}, Bs={self.magnetic_field}, M={self.star_mass}, Mdot={self.star_mass_loss}, D={self.distance} \n"


exoplanets = []
for i in range(150):
    a = 4 * rng.random() - 2  # logAU
    Bs = rng.normal(10, 4)  # gauss
    # Bs = rng.random() * 20  # gauss
    M = rng.normal(0.8, 0.3)  # solar masses
    Mdot = rng.normal(68, 20)  # e-15 solar masses / yr
    D = rng.normal(1000, 400)  # light years
    D = D * 9.46 * 10 ** 15  # conversion to meters
    exo = Exoplanet(a, Bs, M, Mdot, D)
    exoplanets.append(exo)

print(exoplanets)

Bj = 14  # gauss
aj = 5.204  # AU
Rj = 1  # Rj
Rmpj = 60  # Rj
M_sun = 1  # M_sun
M_sun_loss = 68  # e-15 solar masses
nuJ = 24  # MHz

v_sw = 400  # km/s
B_sw = 10 ** (-8)  # teslas

delta = 10.5

frequencies = []
intensities = []
distances = []
semis = []

for exo in exoplanets:
    a = exo.semi_major_axis
    B = exo.magnetic_field
    M = exo.star_mass
    Mdot = exo.star_mass_loss
    D = exo.distance

    semis.append(a)
    a = 10 ** a

    distances.append(D / (9.46 * 10 ** 15))

    nu = nuJ * B / Bj
    frequencies.append(nu)

    Rmp = 56 * Rj * ((B / Bj) / (aj / a)) ** (1 / 3)
    # Lj = 7.4 * 10**22 * Bj**(2/3) * a**(2/3) * v_sw * B_sw**2
    Lj = 2.1 * 10 ** 11

    L = (Rmp / Rmpj) ** 2 * (aj / a) ** 2 * (M / M_sun) * (Mdot / M_sun_loss) * Lj  # Ashtari 2022
    # print(L)

    # I = delta*L / (4*np.pi*D)**2 / nu * 10**(26)
    I = delta * L / (4 * np.pi * nu * D ** 2) * 10 ** 26  # Ashtari 2022
    intensities.append(I)

frequencies = np.array(frequencies)

distances = np.array(distances)
distances = np.reciprocal(distances)
distances *= 10 ** 4

semis = np.array(semis)

df = pd.DataFrame({"x": frequencies,
                   "y": intensities,
                   "d": distances,
                   "s": semis})

fig, ax = plt.subplots()

im = ax.scatter(df.x, df.y, c=df.s, s=df.d, cmap="jet_r")
plt.axvline(x=10, color="black", linestyle="dashed")

fig.colorbar(im, ax=ax, label="Distance to Host Star ($\log_{10}{\mathrm{(AU)}}$)")

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(left=6)

ax.axvspan(0, 10, alpha=0.2, color="teal")

ax.set_title("Frequency and Intensity of CMI Emissions of a Synthetic Exoplanet Sample")
ax.set_xlabel("Emission Frequency (MHz)")
ax.set_ylabel("Radio Brightness (Jy)")
plt.show()

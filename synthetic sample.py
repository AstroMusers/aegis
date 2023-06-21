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
for i in range(100):
    a = rng.normal(0.5, 0.1)  # AU
    Bs = rng.normal(15, 5)  # gauss
    # Bs = rng.random() * 20  # gauss
    M = rng.normal(0.8, 0.3)  # solar masses
    Mdot = rng.normal(68, 10)  # e-15 solar masses
    D = rng.normal(1000, 400)  # light years
    D = D* 9.46 * 10**15
    exo = Exoplanet(a, Bs, M, Mdot, D)
    exoplanets.append(exo)

print(exoplanets)

Bj = 14  # gauss
aj = 5.204  # AU
Rj = 1  # Rj
Rmpj = 60  # Rj
M_sun = 1  # M_sun
M_sun_loss = 68  # e-15 solar masses
nuJ = 24 #MHz

v_sw = 400  # km/s
B_sw = 10 ** (-8)  # teslas

delta = 10.5

frequencies = []
intensities = []

for exo in exoplanets:

    a = exo.semi_major_axis
    B = exo.magnetic_field
    M = exo.star_mass
    Mdot = exo.star_mass_loss
    D = exo.distance

    nu = nuJ * B / Bj
    frequencies.append(nu)

    Rmp = 56*Rj* ( (B/Bj) / (aj/a) )**(1/3)
    # Lj = 7.4 * 10**22 * Bj**(2/3) * a**(2/3) * v_sw * B_sw**2
    Lj = 2.1 * 10**11

    L = (Rmp / Rmpj)**2 * (aj / a)**2 * (M/M_sun) * (Mdot/M_sun_loss) * Lj
    # print(L)

    # I = delta*L / (4*np.pi*D)**2 / nu * 10**(26)
    I = delta*L / (4*np.pi*nu*D)**2 * 10**(26)
    intensities.append(I)

frequencies = np.array(frequencies)
frequencies *= 10**6

df = pd.DataFrame({"x": frequencies,
                   "y": intensities})

plt.scatter(df.x, df.y)
plt.yscale("log")
plt.xscale("log")
plt.ylabel("Radio Intensity (Jy)")
plt.xlabel("Emission Frequency (Hz)")
plt.title("Frequency and Intensity of CMI Emissions of a Synthetic Exoplanet Sample")
plt.show()

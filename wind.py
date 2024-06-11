import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import csv
import pandas as pd
from radio_module import *
import smplotlib


# Read CSV using pandas
filename = "NASA2903.csv"
data = pd.read_csv(filename, comment="#")

# Extract columns into separate arrays
names = data['pl_name'].to_numpy()
semis = data['pl_orbsmax'].to_numpy()
masses = data['st_mass'].to_numpy()
ages = data['st_age'].to_numpy()
spects = data["st_spectype"].to_numpy()
specs = np.array([s[0] for s in spects])


spectral_color = {"M": "xkcd:dark red",
                  "K": "xkcd:tomato red",
                  "G": "xkcd:mango",
                  "F": "xkcd:sunflower yellow"
                  }


k_b = 1.380649 * 10 ** (-29)  # Boltzmann Constant in kg km2 / s2 K
m_p = 1.67262192 * 10 ** (-27)  # Proton mass in kg
G = 8.87129 * 10**2  # Gravitational constant in AU (M_sun)-1 (km/s)2


def non_decreasing(x):
    dx = np.diff(x)
    return np.all(dx >= 0)


def v_at_1_au(t):
    """
    Taken from Lynch2018, citing Newkirk1980
    :param t: Age of the star in yr
    :return: 1 AU speed of the stellar wind in km/s
    """
    v0 = 3.971 * 10**3
    tau = (2.56 * 10 ** (-2))
    return v0 * (1 + t / tau) ** (-0.43)


def sound_speed(T):
    """
    Parker Model
    :param T: Coronal Temperature in K
    :return: sound speed in km/s
    """
    return np.sqrt(k_b * T / m_p)


def critical_radius(M, T):
    """
    Parker Model
    :param M: Stellar mss in M_sun
    :param T: Coronal Temperature in K
    :return: critical radius in AU
    """
    return m_p * G * M / (4 * k_b * T)


targets = np.vectorize(v_at_1_au)(ages)

temperatures = []

for i in range(len(ages)):
    v = targets[i]
    M = masses[i]
    t = ages[i]

    def func(T):
        T = 10**T
        return (v**2 * m_p / (k_b * T)) - np.log((v**2 * m_p / (k_b * T))) - 4 * np.log((4 * k_b * T * 1) / (m_p * G * M)) - (m_p * G * M) / (4 * k_b * T * 1) + 3

    guess = np.array([6.74])
    soln = fsolve(func, guess)
    actual = 10**soln[0]

    temperatures.append(actual)

radial_ranges = []
wind_speeds = []

for j in range(len(ages)):
    M = masses[j]
    t = ages[j]
    T = temperatures[j]
    a = semis[j]

    r_c = critical_radius(M, T)

    def fn(v, r):
        v = 10**v
        # r = a
        return (v**2 * m_p / (k_b * T)) - np.log((v**2 * m_p / (k_b * T))) - 4 * np.log((4 * k_b * T * r) / (m_p * G * M)) - 4 * (m_p * G * M) / (4 * k_b * T * r) + 3

    cond = True

    def v_for_r(r):

        fast_bois_1512 = [792]
        fast_bois_1912 = [1940, 3643]
        fast_bois_3112 = [1164]
        fast_bois_2903 = [13, 243, 440, 1350, 1356]
        if r < r_c:
            guess = 0
        else:
            guess = 3
            if j in fast_bois_2903:
                guess = 3.5

        guess = np.array([guess])

        def func(v):
            return fn(v, r=r)

        soln = fsolve(func, guess)
        return 10**soln[0]

    # r_values = np.linspace(0.1*r_c, 200*r_c, 100)
    r_values = np.logspace(np.log10(0.1*r_c), np.log10(200*r_c), 100)
    v_values = [v_for_r(r) for r in r_values]

    radial_ranges.append(r_values)
    wind_speeds.append(v_values)


def profile_plot(rs, ws, save=False):
    plt.rcParams['figure.figsize'] = [5, 4]
    plt.rcParams["font.size"] = 13
    fig, ax = plt.subplots()

    plotted_spectral_types = set()

    for i in range(len(rs)):
        spectral_type = specs[i]
        crit_index = np.argmin(abs(rs[i] - rs[i][0]*10))
        ax.plot(rs[i][crit_index], ws[i][crit_index], "k*", markersize=3)
        if i == 0:
            ax.plot([], [], "k*", label="$r_c$", markersize=10)
        color = spectral_color.get(spectral_type, "black")
        ax.plot(rs[i], ws[i], color=color, alpha=0.5)
        if spectral_type not in plotted_spectral_types:
            plt.plot([], [], label=spectral_type, color=color)
            plotted_spectral_types.add(spectral_type)
        # ax.scatter(rs[i], ws[i], color="k", alpha=0.05)
        ax.set(yscale="log", xscale="log")
        ax.set_ylim(10, 5000)
        ax.set_xlim(3e-5, 50)
        ax.set(xlabel="Distance from Host Star (AU)",
               ylabel="Parker Wind Speed (km/h)")
    # ax.minorticks_on()
    # retro_noir(ax)
    plt.legend()
    plt.tight_layout()
    plt.show()
    if save:
        plt.savefig("wind_profiles.pdf")


profile_plot(radial_ranges, wind_speeds)

# plt.show()
problem = False
for i in range(len(wind_speeds)):
    profile = wind_speeds[i]
    if not non_decreasing(profile[5:]):
        problem = True
        print(f"Problem with Profile {i}!")

if not problem:
    print("All good, fantastic work!")
    speeds_at_distance = []
    for i in range(len(wind_speeds)):
        a = semis[i]
        v = v_for_r(a)
        speeds_at_distance.append(v)

temperatures = np.round(np.array(temperatures))
speeds_at_distance = np.round(np.array(speeds_at_distance))
results = np.vstack((names, temperatures, speeds_at_distance))
results = results.T

column_headers = ["Name", "Corona Temp. (K)", "Wind Speed (km/s)"]
windfile = "wind_info" + filename[4:-4] + ".txt"
with open(windfile, mode="w") as file:
    writer = csv.writer(file, delimiter="\t")
    writer.writerow(column_headers)
    writer.writerows(results)

print(f"Results saved to {windfile}")

# plt.show()
# Code takes 19 seconds to run end to end when radial profile linspace has 100 elements

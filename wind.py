import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.optimize import root

semi_column = 6
age_column = 32
st_mass_column = 28

names, semis, masses, ages = np.genfromtxt("NASA0808.csv", usecols=(0, semi_column, st_mass_column, age_column), skip_header=56,
                                    filling_values=0, delimiter=",", unpack=True)

k_b = 1.380649 * 10 ** (-29)  # Boltzmann Constant in kg km2 / s2 K
m_p = 1.67262192 * 10 ** (-27)  # Proton mass in kg
G = 8.87129 * 10**2  # Gravitational constant in AU (M_sun)-1 (km/s)2


def non_decreasing(x):
    dx = np.diff(x)
    return np.all(dx >= 0)


def v_at_1_AU(t):
    """

    :param t: Age of the star in yr
    :return: 1 AU speed of the stellar wind in km/s
    """
    v0 = 3.971 * 10**3
    tau = (2.56 * 10 ** 7)
    return v0 * (1 + t / tau) ** (-0.43)


def sound_speed(T):
    """

    :param T: Coronal Temperature in K
    :return: sound speed in km/s
    """
    return np.sqrt(k_b * T / m_p)


def critical_radius(M, T):
    """

    :param M:
    :param T:
    :return: critical radius in AU
    """
    return m_p * G * M / (4 * k_b * T)


ages *= 10**9
targets = np.vectorize(v_at_1_AU)(ages)
print(ages)
print(masses)
print(targets)

# plt.hist(targets)
# plt.show()

trg_lst = targets.tolist()
age_lst = ages.tolist()
mas_lst = masses.tolist()

temperatures = []

for i in range(len(trg_lst)):
    v = trg_lst[i]
    M = masses[i]
    t = ages[i] * 10**(-9)
    print(v, M, t)

    def func(T):
        T = 10**T
        return (v**2 * m_p / (k_b * T)) - np.log((v**2 * m_p / (k_b * T))) - 4 * np.log((4 * k_b * T * 1) / (m_p * G * M)) - (m_p * G * M) / (4 * k_b * T * 1) +3

    guess = 6.74
    print(guess)
    soln = fsolve(func, guess)
    actual = 10**soln[0]

    temperatures.append(actual)

print(temperatures)

radial_ranges = []
wind_speeds = []

for j in range(len(age_lst)):
    M = masses[j]
    t = ages[j]
    T = temperatures[j]
    a = semis[j]

    r_c = critical_radius(M, T)

    def fn(v, r):
        v = 10**v
        # r = a
        return (v**2 * m_p / (k_b * T)) - np.log((v**2 * m_p / (k_b * T))) - 4 * np.log((4 * k_b * T * r) / (m_p * G * M)) - (m_p * G * M) / (4 * k_b * T * r) +3

    cond = True
    def v_for_r(r):

        guess = 3
        if j == 792:
            guess = 3.5

        def func(v):
            return fn(v, r=r)

        soln = fsolve(func, guess)
        return 10**soln[0]

    # upper = max(a + 1, 200*r_c)
    r_values = np.linspace(0.1*r_c, 200*r_c, 100)
    v_values = [v_for_r(r) for r in r_values]

    radial_ranges.append(r_values)
    wind_speeds.append(v_values)

    # print(r_values)
    # print(v_values)

    plt.plot(r_values, v_values, alpha=0.2)

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

print(speeds_at_distance)


# Code takes 19 seconds to run end to end when radial profile linspace has 100 elements

plt.show()
# print(wind_speeds)


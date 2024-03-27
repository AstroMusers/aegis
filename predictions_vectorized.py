import numpy as np
import pandas as pd
import math
from tabulate import tabulate
from adjustText import adjust_text
from radio_module import *
from rotation_script import *

rng = np.random.default_rng()

# Retrieve Data
filename = "NASA3112.csv"
df = pd.read_csv(filename, comment="#")
windfile = "wind_info" + filename[4:-4] + ".txt"


def calculate_mass(row):
    if row['pl_bmassprov'] in ["Mass", "Msin(i)/sin(i)"]:
        M_i = row['pl_bmassj']
    else:
        M_i = row['pl_bmassj'] * 1.15  # Expected Value of the mass based on projected mass
    M_ss = ((row['pl_bmassjerr1'] - row['pl_bmassjerr2']) / 2)
    return M_i, M_ss


def reject_sample(values, std):
    sampled_values = rng.normal(values, std)
    while np.any(sampled_values <= 0):
        mask = sampled_values <= 0
        sampled_values[mask] = rng.normal(values[mask], std[mask])
    return sampled_values


wind_temperatures, wind_speeds = np.genfromtxt(windfile, usecols=(1, 2), skip_header=1, delimiter="\t", unpack=True)

df["orbper_std"] = ((df["pl_orbpererr1"] - df["pl_orbpererr2"]) / 2).fillna(df["pl_orbper"] / 5)
df["orbsmax_std"] = ((df["pl_orbsmaxerr1"] - df["pl_orbsmaxerr2"]) / 2).fillna(df["pl_orbsmax"] / 5)
df["pl_rad_std"] = ((df["pl_radjerr1"] - df["pl_radjerr2"]) / 2).fillna(df["pl_radj"] / 5)
df[['pl_mass', 'pl_mass_std']] = df.apply(calculate_mass, axis=1, result_type='expand')
df["pl_mass_std"] = df["pl_mass_std"].fillna(df["pl_mass"] / 5)
df["pl_dens_std"] = ((df["pl_denserr1"] - df["pl_denserr2"]) / 2).fillna(df["pl_dens"] / 5)
df["st_mass_std"] = ((df["st_masserr1"] - df["st_masserr2"]) / 2).fillna(df["st_mass"] / 5)
df["st_age_std"] = ((df["st_ageerr1"] - df["st_ageerr2"]) / 2).fillna(df["st_age"] / 5)
df["st_rad_std"] = ((df["st_raderr1"] - df["st_raderr2"]) / 2).fillna(df["st_rad"] / 5)
# df["T_wind"] = wind_temperatures
# df["v_wind"] = wind_speeds
T_wind = wind_temperatures
v_wind = wind_speeds
df["d_ly"] = df["sy_dist"] * 3.261561

rotation = np.vectorize(rotation)
convective_radius = np.vectorize(convective_radius)
Rm = np.vectorize(Rm)

for k in range(1000):
    T = reject_sample(df["pl_orbper"], df["orbper_std"])
    a = reject_sample(df["pl_orbsmax"], df["orbsmax_std"])
    R = reject_sample(df["pl_radj"], df["pl_rad_std"])
    M = reject_sample(df["pl_mass"], df["pl_mass_std"])
    p = reject_sample(df["pl_dens"], df["pl_dens_std"])
    M_s = reject_sample(df["st_mass"], df["st_mass_std"])
    t = reject_sample(df["st_age"], df["st_age_std"])
    Rs = reject_sample(df["st_rad"], df["st_rad_std"])

    b0_exponent = rng.normal(-0.655, 0.045)  # Vidotto 2014
    flux_exponent = rng.normal(-1.74, 0.34)  # Ayres 1997 in Lynch 2018
    loss_exponent = rng.normal(0.79, 0.17)  # Alvarado-Gomez 2016 in Lynch 2018

    L = moment_sampler()

    highS_Mdot = t ** (-1.23) * 10 ** 3
    lowS_Mdot = t ** (-0.9) * 10 ** 3

    Mdot = mass_loss(t, flux_exponent, loss_exponent)

    sigma = 1  # Jupiter conductivity

    p_c = density(p)
    w_p = rotation(T, a, L, M, R)

    r_c = convective_radius(M, p_c, R)
    mu = magnetic_moment(p_c, w_p, r_c, sigma)
    B = magnetic_field(mu, R)

    D = (df["d_ly"] * 9.46 * 10 ** 15).values

    B_perp = imf_perp_complete(M_s, a, Rs, t, v_wind, b0_exponent)
    v_k = keplerian(M_s, a)
    veff = v_eff(v_wind, v_k)
    B_star = imf_complete(M_s, a, v_wind, t, Rs, b0_exponent)
    n = number_density(Mdot, veff, a)
    R_m = Rm(B, R, n, T_wind, v_wind, B_star)

    if k == 999:
        print("wtf how fast")



from pathlib import Path
from radio_module import *
from rotation_script import *
import pandas as pd
import numpy as np
from repository_paths import get_repository_path
from tdpy.verbosity import print

DATA_DIR = get_repository_path()
filename = str(DATA_DIR / "NASA3004.csv")
print(f"Reading from {filename}...")
df = pd.read_csv(filename, comment="#")
windfile = DATA_DIR / f"stellar_wind_predictions_{Path(filename).stem[4:]}.tsv"

print(f"Reading from {windfile}...")
wind_temperatures, wind_speeds = np.genfromtxt(windfile, usecols=(1, 2), skip_header=1, delimiter="\t", unpack=True)


def default_calculate(name, burst=10, filename=filename):
    print(f"Reading from {filename}...")
    df = pd.read_csv(filename, comment="#")
    ind = df.index[df["pl_name"] == name].tolist()[0]
    # print(ind)
    T, a, R, p, M_s, t, Rs, = df.loc[ind, ["pl_orbper", "pl_orbsmax", "pl_radj", "pl_dens", "st_mass", "st_age", "st_rad"]]
    print(T, a, R)
    T_wind, v_wind = wind_temperatures[ind], wind_speeds[ind]
    if (df.loc[ind, "pl_bmassprov"] == "Mass") or (df.loc[ind, "pl_bmassprov"] == "Msin(i)/sin(i)"):
        M = df.loc[ind, "pl_bmassj"]
    else:
        M = df.loc[ind, "pl_bmassj"] * 1.15  # Expected Value of the mass based on projected mass
    b0_exponent = -0.655
    flux_exponent = -1.74
    loss_exponent = 0.79
    d = df.loc[ind, "sy_dist"] * 3.261561
    L = moment_sampler()
    Mdot = mass_loss(t, flux_exponent, loss_exponent)
    sigma = 1  # Jupiter conductivity
    p_c = density(p)
    w_p = rotation(T, a, L, M, R)
    r_c = convective_radius(M, p_c, R)
    mu = magnetic_moment(p_c, w_p, r_c, sigma)
    B = magnetic_field(mu, R)
    D = d * 9.46 * 10 ** 15  # conversion to meters
    B_perp = imf_perp_complete(M_s, a, Rs, t, v_wind, b0_exponent)
    v_k = keplerian(M_s, a)
    veff = v_eff(v_wind, v_k)
    B_star = imf_complete(M_s, a, v_wind, t, Rs, b0_exponent)
    n = number_density(Mdot, veff, a)
    R_m = Rm(B, R, n, T_wind, v_wind, B_star)
    n_p = 8.98 * np.sqrt(n) * 10 ** (-3)
    nu = max_freq(B)
    n_p /= 10 ** 6
    # I = complete(B, a, M_s, Mdot, D)
    P_in_mag = P_input_mag(B_perp, veff * 10 ** 3, R_m * 7 * 10 ** 8, n)
    P_in_kin = P_input_kin(B_perp, veff * 10 ** 3, R_m * 7 * 10 ** 8, n)
    P_rad_both = radio_power(P_in_mag, P_in_kin, nu, D, both=True) * burst
    print(f"v_max = {nu/1e6}, I={P_rad_both}")
    return float(nu/1e6), float(P_rad_both)


if __name__ == "__main__":
    default_calculate("tau Boo b")
    calculator = np.vectorize(default_calculate)

    predictions_path = DATA_DIR / "Output Tables" / "all.csv"
    print(f"Reading from {predictions_path}...")
    MCdf = pd.read_csv(predictions_path)
    meta_sub = df[df["pl_orbsmax"] < 0.075][["pl_name", "pl_orbsmax", "pl_orbsmaxerr1"]]
    meta_sub.rename(columns={"pl_name": "Name"}, inplace=True)
    merged = pd.merge(MCdf, meta_sub, on="Name")
    merged["Freq"], merged["Flux"] = zip(*merged["Name"].apply(calculator))
    merged["freq_ratio"] = merged["Max. Frequency [MHz]"] / merged["Freq"]
    merged["flux_ratio"] = merged["Max. Flux Density [Jy]"] / merged["Flux"]
    freq_sorted = merged.sort_values(by="freq_ratio", ascending=False)
    flux_sorted = merged.sort_values(by="flux_ratio", ascending=False)
    print(len(merged["freq_ratio"][(merged["freq_ratio"] < 0.667) | (merged["freq_ratio"] > 1.5)]))
    print(len(merged["flux_ratio"][(merged["flux_ratio"] < 0.667) | (merged["flux_ratio"] > 1.5)]))

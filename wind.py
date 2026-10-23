from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import pandas as pd
from radio_module import *
from repository_paths import get_repository_path

DATA_DIR = get_repository_path()


# Read CSV using pandas
filename = str(DATA_DIR / "NASA3004.csv")
data = pd.read_csv(filename, comment="#")

# Extract columns into separate arrays
names = data['pl_name'].to_numpy()
semis = data['pl_orbsmax'].to_numpy()
masses = data['st_mass'].to_numpy()
ages = data['st_age'].to_numpy()
spects = data["st_spectype"].fillna("NA").to_numpy()
specs = np.array([s[0] for s in spects])


spectral_color = {"M": "xkcd:dark red",
                  "K": "xkcd:tomato red",
                  "G": "xkcd:mango",
                  "F": "xkcd:sunflower yellow"
                  }


def profile_plot(dfs, save=False, show=False):

    fig, axs = plt.subplots(2, 2, figsize=(6, 6), sharex=True, sharey=True)

    plotted_spectral_types = set()

    for k, (type, df) in enumerate(dfs.items()):
        ax = axs[k // 2, k % 2]
        for i in range(len(df.rs)):
            crit_index = np.argmin(abs(df.rs[i] - df.rs[i][0]*10))
            ax.plot(df.rs[i][crit_index], df.ws[i][crit_index], "k*", markersize=3, alpha=0.5)
            if i == 0 and k == 0:
                ax.plot([], [], "k*", label="$r_c$", markersize=10)
            color = spectral_color.get(type, "black")
            ax.plot(df.rs[i], df.ws[i], color=color, alpha=0.1)
            if type not in plotted_spectral_types:
                plt.plot([], [], label=type, color=color)
                plotted_spectral_types.add(type)
        ax.set(yscale="log", xscale="log")
        ax.grid(True, alpha=0.5)
        ax.set_ylim(10, 5000)
        ax.set_xlim(3e-5, 50)
        retro_noir(ax)

    fig.supylabel("Parker Wind Speed (km/h)")
    fig.supxlabel("Distance From Host Star (AU)")
    fig.suptitle("\n")
    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    fig.legend(lines, labels, loc="upper center", bbox_to_anchor=(0.5, 1), ncol=5, frameon=True, shadow=True)
    plt.subplots_adjust(top=0.85)  # Increase top margin to make space for the legend
    plt.tight_layout()  # Adjust rect parameter to leave space for the legend

    if save:
        plt.savefig(DATA_DIR / "winds.pdf")

    if show:
        plt.show()
    else:
        plt.close(fig)


def run_wind_simulation(save=True, show=False):
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
    speed_functions = []

    for j in range(len(ages)):
        M = masses[j]
        t = ages[j]
        T = temperatures[j]
        a = semis[j]

        r_c = critical_radius(M, T)

        def fn(v, r, T, M):
            v = 10**v
            return (v**2 * m_p / (k_b * T)) - np.log((v**2 * m_p / (k_b * T))) - 4 * np.log((4 * k_b * T * r) / (m_p * G * M)) - 4 * (m_p * G * M) / (4 * k_b * T * r) + 3

        def v_for_r(r, T=T, M=M):
            fast_bois_3004 = [382, 1362]

            if r < r_c:
                guess = 0
            else:
                guess = 3
                if j in fast_bois_3004:
                    guess = 3.5

            guess = np.array([guess])

            def func(v):
                return fn(v, r=r, T=T, M=M)

            soln = fsolve(func, guess)
            return 10**soln[0]

        speed_functions.append(v_for_r)

        r_values = np.logspace(np.log10(0.1*r_c), np.log10(200*r_c), 100)
        v_values = [v_for_r(r) for r in r_values]

        radial_ranges.append(r_values)
        wind_speeds.append(v_values)

    data1 = {"names": names,
             "rs": radial_ranges,
             "ws": wind_speeds,
             "specs": specs}
    df = pd.DataFrame(data1)
    df["star_names"] = df["names"].str[:-1]
    df = df.drop_duplicates(subset=["star_names"])
    df = df[df["specs"] != "N"]
    df = df[df["specs"] != "A"]

    grouped = df.groupby("specs")
    dfs = {category: group.reset_index(drop=True) for category, group in grouped}

    profile_plot(dfs, save=save, show=show)

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
            speed_function = speed_functions[i]
            v = speed_function(a)
            speeds_at_distance.append(v)

    temperatures = np.round(np.array(temperatures))
    speeds_at_distance = np.round(np.array(speeds_at_distance))
    results = np.vstack((names, temperatures, speeds_at_distance))
    results = results.T

    column_headers = ["Name", "Corona Temp. (K)", "Wind Speed (km/s)"]
    windfile = str(
        DATA_DIR / f"stellar_wind_predictions_{Path(filename).stem[4:]}.tsv"
    )
    with open(windfile, mode="w") as file:
        writer = csv.writer(file, delimiter="\t")
        writer.writerow(column_headers)
        writer.writerows(results)

    print(f"Results saved to {windfile}")

    index = int(np.where(names == "tau Boo b")[0])
    np.savez(DATA_DIR / "taub_wind", ranges=radial_ranges[index], speeds=wind_speeds[index])


if __name__ == "__main__":
    run_wind_simulation(save=True, show=True)

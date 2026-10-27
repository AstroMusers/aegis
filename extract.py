from pathlib import Path
import numpy as np
import pandas as pd
from radio_module import *

from coordinate_formatting import build_radio_catalog_table, three_format, two_format
from repository_paths import get_repository_path
from tdpy.verbosity import print

DATA_DIR = get_repository_path()


def extract_observable_targets(save=True):
    filename = DATA_DIR / "NASA3004.csv"
    print(f"Reading from {filename}...")
    df = pd.read_csv(filename, comment="#")

    detectable_path = DATA_DIR / "detectables_both.csv"
    print(f"Reading from {detectable_path}...")
    df2 = pd.read_csv(detectable_path)
    final = build_radio_catalog_table(df, df2, "Freq(MHz)", "Flux(mJy)")

    final["Mass (MJ)"] = final["Mass (MJ)"].apply(three_format)
    final["a (AU)"] = final["a (AU)"].apply(three_format)

    final["Radius (RJ)"] = final["Radius (RJ)"].apply(two_format)
    final["d (pc)"] = final["d (pc)"].apply(two_format)
    final["t (Gyr)"] = final["t (Gyr)"].apply(two_format)
    final["vmax"] = final["vmax"].apply(two_format)

    if save:
        output_path = DATA_DIR / "obs_table.csv"
        print(f"Writing to {output_path}...")
        final.to_csv(output_path, index=False)

    final["L"] = final["Phi (mJy)"] * final["d (pc)"].apply(float) **2
    more_info = final.sort_values(by="L", ascending=False).copy()
    more_info["R_earth"] = more_info["Radius (RJ)"].apply(float) * 11.209
    more_info["M_earth"] = more_info["Mass (MJ)"].apply(float) * 318
    return final, more_info


if __name__ == "__main__":
    extract_observable_targets(save=True)

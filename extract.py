from pathlib import Path
import numpy as np
import pandas as pd
from radio_module import *

from coordinate_formatting import format_coords, three_format, two_format
from repository_paths import get_repository_path

DATA_DIR = get_repository_path()


def extract_observable_targets(save=True):
    filename = DATA_DIR / "NASA3004.csv"
    df = pd.read_csv(filename, comment="#")

    df2 = pd.read_csv(DATA_DIR / "detectables_both.csv")
    detect_data = df[df["pl_name"].isin(df2["Name"])]
    enough_data = detect_data[["pl_name", "pl_bmassj", "pl_radj", "pl_orbsmax", "sy_dist", "st_age", "ra", "dec"]].copy()
    enough_data.rename(columns={"pl_name": "Name"}, inplace=True)

    merged = pd.merge(enough_data, df2, on="Name")

    merged[['ra', 'dec']] = merged.apply(format_coords, axis=1)

    merged["vmax"] = merged["Freq(MHz)"]
    merged.columns = ["Name", "Mass (MJ)", "Radius (RJ)", "a (AU)", "d (pc)", "t (Gyr)", "RA (J2000)", "DEC (J2000)", "nupeak (MHz)", "Phi (mJy)", "vmax"]
    final = merged[["Name", "RA (J2000)", "DEC (J2000)", "Mass (MJ)", "Radius (RJ)", "a (AU)", "d (pc)", "t (Gyr)", "vmax", "Phi (mJy)"]]
    final = final.sort_values(by="Phi (mJy)", ascending=False)

    final["Mass (MJ)"] = final["Mass (MJ)"].apply(three_format)
    final["a (AU)"] = final["a (AU)"].apply(three_format)

    final["Radius (RJ)"] = final["Radius (RJ)"].apply(two_format)
    final["d (pc)"] = final["d (pc)"].apply(two_format)
    final["t (Gyr)"] = final["t (Gyr)"].apply(two_format)
    final["vmax"] = final["vmax"].apply(two_format)

    if save:
        final.to_csv(DATA_DIR / "obs_table.csv", index=False)

    final["L"] = final["Phi (mJy)"] * final["d (pc)"].apply(float) **2
    more_info = final.sort_values(by="L", ascending=False).copy()
    more_info["R_earth"] = more_info["Radius (RJ)"].apply(float) * 11.209
    more_info["M_earth"] = more_info["Mass (MJ)"].apply(float) * 318
    return final, more_info


if __name__ == "__main__":
    extract_observable_targets(save=True)

from pathlib import Path
import pandas as pd

from coordinate_formatting import format_coords, one_format, three_format, two_format, zero_format
from repository_paths import get_repository_path
from tdpy.verbosity import print

DATA_DIR = get_repository_path()


def append_predictions(save=True):
    filename = DATA_DIR / "NASA3004.csv"
    print(f"Reading from {filename}...")
    df = pd.read_csv(filename, comment="#")

    predictions_path = DATA_DIR / "Output Tables" / "all.csv"
    print(f"Reading from {predictions_path}...")
    df2 = pd.read_csv(predictions_path)
    print(len(df2[0.1 < df2["Max. Frequency [MHz]"]]))
    df2 = df2[(0.1 < df2["Max. Frequency [MHz]"]) & (df2["Max. Frequency [MHz]"] < 10)]
    if save:
        filtered_path = DATA_DIR / "Output Tables" / "0.1-10MHz.csv"
        print(f"Writing to {filtered_path}...")
        df2.to_csv(filtered_path, index=False)
    df2 = df2[:5]
    df2 = df2[["Name", "Max. Frequency [MHz]", "Max. Flux Density [Jy]"]]

    detect_data = df[df["pl_name"].isin(df2["Name"])]
    enough_data = detect_data[["pl_name", "pl_bmassj", "pl_radj", "pl_orbsmax", "sy_dist", "st_age", "ra", "dec"]].copy()
    enough_data.rename(columns={"pl_name": "Name"}, inplace=True)

    merged = pd.merge(enough_data, df2, on="Name")

    merged[['ra', 'dec']] = merged.apply(format_coords, axis=1)
    merged["vmax"] = merged["Max. Frequency [MHz]"]
    merged.columns = ["Name", "Mass (MJ)", "Radius (RJ)", "a (AU)", "d (pc)", "t (Gyr)", "RA (J2000)", "DEC (J2000)", "nupeak (MHz)", "Phi (mJy)", "vmax"]
    final = merged[["Name", "RA (J2000)", "DEC (J2000)", "Mass (MJ)", "Radius (RJ)", "a (AU)", "d (pc)", "t (Gyr)", "vmax", "Phi (mJy)"]]
    final = final.sort_values(by="Phi (mJy)", ascending=False)

    final["Mass (MJ)"] = final["Mass (MJ)"].apply(three_format)
    final["a (AU)"] = (final["a (AU)"]*1e2).apply(one_format)

    final["Radius (RJ)"] = final["Radius (RJ)"].apply(two_format)
    final["d (pc)"] = final["d (pc)"].apply(zero_format)
    final["t (Gyr)"] = final["t (Gyr)"].apply(two_format)
    final["vmax"] = final["vmax"].apply(two_format)
    final["Phi (Jy)"] = (final["Phi (mJy)"]).apply(two_format)

    if save:
        output_path = DATA_DIR / "append_table.csv"
        print(f"Writing to {output_path}...")
        final.to_csv(output_path, columns=["Name", "RA (J2000)", "DEC (J2000)", "vmax", "Phi (Jy)"], index=False)
    return final


if __name__ == "__main__":
    append_predictions(save=True)

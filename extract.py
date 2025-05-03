import numpy as np
import pandas as pd
from radio_module import *
from astropy.coordinates import SkyCoord
import astropy.units as u

filename = "NASA3004.csv"
df = pd.read_csv(filename, comment="#")

# det, freq, flux = np.genfromtxt("detectables_both.txt", dtype=str, delimiter="  ", skip_header=1, skip_footer=1, usecols=(0,1,2), unpack=True)

df2 = pd.read_csv("detectables_both.csv")
detect_data = df[df["pl_name"].isin(df2["Name"])]
enough_data = detect_data[["pl_name", "pl_bmassj", "pl_radj", "pl_orbsmax", "sy_dist", "st_age", "ra", "dec"]]
enough_data.rename(columns={"pl_name": "Name"}, inplace=True)

merged = pd.merge(enough_data, df2, on="Name")
# merged["System Name"] = merged["Name"]


def two_format(x):
    return '{:.2f}'.format(x)


def three_format(x):
    return '{:.3f}'.format(x)


# Function to format RA/DEC
def format_coords(row):
    coord = SkyCoord(ra=row['ra']*u.deg, dec=row['dec']*u.deg)
    ra_str = coord.ra.to_string(unit=u.hour, sep=':', precision=0, pad=True)
    dec_str = coord.dec.to_string(unit=u.deg, sep=':', precision=0, alwayssign=True, pad=True)
    return pd.Series({'RA': ra_str, 'DEC': dec_str})

merged[['ra', 'dec']] = merged.apply(format_coords, axis=1)

merged["vmax"] = merged["Freq(MHz)"]
merged.columns = ["Name", "Mass (MJ)", "Radius (RJ)", "a (AU)", "d (pc)", "t (Gyr)", "RA (J2000)", "DEC (J2000)", "nupeak (MHz)", "Phi (mJy)", "vmax"]
final = merged[["Name", "RA (J2000)", "DEC (J2000)", "Mass (MJ)", "Radius (RJ)", "a (AU)", "d (pc)", "t (Gyr)", "vmax", "Phi (mJy)"]]
final = final.sort_values(by="Phi (mJy)", ascending=False)

# final["RA (J2000)"] = final["RA (J2000)"].apply(format_ra)
# final["DEC (J2000)"] = final["DEC (J2000)"].apply(format_ra)

final["Mass (MJ)"] = final["Mass (MJ)"].apply(three_format)
final["a (AU)"] = final["a (AU)"].apply(three_format)

final["Radius (RJ)"] = final["Radius (RJ)"].apply(two_format)
final["d (pc)"] = final["d (pc)"].apply(two_format)
final["t (Gyr)"] = final["t (Gyr)"].apply(two_format)
final["vmax"] = final["vmax"].apply(two_format)

final.to_csv("obs_table.csv", index=False)

final["L"] = final["Phi (mJy)"] * final["d (pc)"].apply(float) **2
more_info = final.sort_values(by="L")
more_info = final.sort_values(by="L", ascending=False)
more_info["R_earth"] = more_info["Radius (RJ)"].apply(float) * 11.209
more_info["M_earth"] = more_info["Mass (MJ)"].apply(float) * 318

# enough_data["freq"], detect_data["flux"] = freq, flux
# sorted = enough_data.sort_values(by="pl_name")
import numpy as np
import pandas as pd
from radio_module import *
from astroquery.simbad import Simbad

filename = "NASA2903.csv"
df = pd.read_csv(filename, comment="#")

# det, freq, flux = np.genfromtxt("detectables_both.txt", dtype=str, delimiter="  ", skip_header=1, skip_footer=1, usecols=(0,1,2), unpack=True)

df2 = pd.read_csv("detectables_both.csv")
detect_data = df[df["pl_name"].isin(df2["Name"])]
enough_data = detect_data[["pl_name", "pl_bmassj", "pl_radj", "pl_orbsmax", "sy_dist", "st_age"]]
enough_data.rename(columns={"pl_name": "Name"}, inplace=True)

merged = pd.merge(enough_data, df2, on="Name")
# merged["System Name"] = merged["Name"]


def two_format(x):
    return '{:.2f}'.format(x)


def three_format(x):
    return '{:.3f}'.format(x)


def query_simbad(object_name):
    result_table = Simbad.query_object(object_name)
    if result_table is not None:
        ra = result_table['RA'][0]
        dec = result_table['DEC'][0]
        return ra, dec
    else:
        return None, None


def format_ra(ra_string):
    hours, minutes, seconds = ra_string.split()
    seconds = round(float(seconds))
    return f"{hours}:{minutes}:{seconds:02d}"


merged['RA'], merged['DEC'] = zip(*(merged["Name"].str[:-1]).apply(query_simbad))
merged["vmax"] = merged["Freq(MHz)"]
merged.columns = ["Name", "Mass (MJ)", "Radius (RJ)", "a (AU)", "d (pc)", "t (Gyr)", "nupeak (MHz)", "Phi (mJy)", "RA (J2000)", "DEC (J2000)", "vmax"]
final = merged[["Name", "RA (J2000)", "DEC (J2000)", "Mass (MJ)", "Radius (RJ)", "a (AU)", "d (pc)", "t (Gyr)", "vmax", "Phi (mJy)"]]
final = final.sort_values(by="Phi (mJy)", ascending=False)

final["RA (J2000)"] = final["RA (J2000)"].apply(format_ra)
final["DEC (J2000)"] = final["DEC (J2000)"].apply(format_ra)

final["Mass (MJ)"] = final["Mass (MJ)"].apply(three_format)
final["a (AU)"] = final["a (AU)"].apply(three_format)

final["Radius (RJ)"] = final["Radius (RJ)"].apply(two_format)
final["d (pc)"] = final["d (pc)"].apply(two_format)
final["t (Gyr)"] = final["t (Gyr)"].apply(two_format)
final["vmax"] = final["vmax"].apply(two_format)

final.to_csv("obs_table.csv", index=False)

final["L"] = final["Phi (mJy)"] * final["d (pc)"].apply(float) **2
zort = final.sort_values(by="L")
zort = final.sort_values(by="L", ascending=False)
zort["R_earth"] = zort["Radius (RJ)"].apply(float) * 11.209

# enough_data["freq"], detect_data["flux"] = freq, flux
# sorted = enough_data.sort_values(by="pl_name")
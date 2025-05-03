import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u

filename = "NASA3004.csv"
df = pd.read_csv(filename, comment="#")

df2 = pd.read_csv("Output Tables/all.csv")
print(len(df2[0.1 < df2["Max. Frequency [MHz]"]]))
df2 = df2[(0.1 < df2["Max. Frequency [MHz]"]) & (df2["Max. Frequency [MHz]"] < 10)]
df2.to_csv("Output Tables/0.1-10MHz.csv", index=False)
df2 = df2[:5]
df2 = df2[["Name", "Max. Frequency [MHz]", "Max. Flux Density [Jy]"]]

detect_data = df[df["pl_name"].isin(df2["Name"])]
enough_data = detect_data[["pl_name", "pl_bmassj", "pl_radj", "pl_orbsmax", "sy_dist", "st_age", "ra", "dec"]]
enough_data.rename(columns={"pl_name": "Name"}, inplace=True)

merged = pd.merge(enough_data, df2, on="Name")

def two_format(x):
    return '{:.2f}'.format(x)


def three_format(x):
    return '{:.3f}'.format(x)


def format_coords(row):
    coord = SkyCoord(ra=row['ra']*u.deg, dec=row['dec']*u.deg)
    ra_str = coord.ra.to_string(unit=u.hour, sep=':', precision=0, pad=True)
    dec_str = coord.dec.to_string(unit=u.deg, sep=':', precision=0, alwayssign=True, pad=True)
    return pd.Series({'RA': ra_str, 'DEC': dec_str})


def one_format(x):
    return '{:.1f}'.format(x)


def zero_format(x):
    return '{:.0f}'.format(x)

merged[['ra', 'dec']] = merged.apply(format_coords, axis=1)
merged["vmax"] = merged["Max. Frequency [MHz]"]
merged.columns = ["Name", "Mass (MJ)", "Radius (RJ)", "a (AU)", "d (pc)", "t (Gyr)", "RA (J2000)", "DEC (J2000)", "nupeak (MHz)", "Phi (mJy)", "vmax"]
final = merged[["Name", "RA (J2000)", "DEC (J2000)", "Mass (MJ)", "Radius (RJ)", "a (AU)", "d (pc)", "t (Gyr)", "vmax", "Phi (mJy)"]]
final = final.sort_values(by="Phi (mJy)", ascending=False)

# final["RA (J2000)"] = final["RA (J2000)"].apply(format_ra)
# final["DEC (J2000)"] = final["DEC (J2000)"].apply(format_ra)

final["Mass (MJ)"] = final["Mass (MJ)"].apply(three_format)
final["a (AU)"] = (final["a (AU)"]*1e2).apply(one_format)

final["Radius (RJ)"] = final["Radius (RJ)"].apply(two_format)
final["d (pc)"] = final["d (pc)"].apply(zero_format)
final["t (Gyr)"] = final["t (Gyr)"].apply(two_format)
final["vmax"] = final["vmax"].apply(two_format)
final["Phi (Jy)"] = (final["Phi (mJy)"]).apply(two_format)

final.to_csv("append_table.csv", columns=["Name", "RA (J2000)", "DEC (J2000)", "vmax", "Phi (Jy)"], index=False)
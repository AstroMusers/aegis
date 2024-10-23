import pandas as pd
from astroquery.simbad import Simbad

filename = "NASA1407.csv"
df = pd.read_csv(filename, comment="#")

df2 = pd.read_csv("Output Tables/all.csv")
print(len(df2[0.1 < df2["Max. Frequency [MHz]"]]))
df2 = df2[(0.1 < df2["Max. Frequency [MHz]"]) & (df2["Max. Frequency [MHz]"] < 10)][:5]
df2 = df2[["Name", "Max. Frequency [MHz]", "Max. Flux Density [Jy]"]]
detect_data = df[df["pl_name"].isin(df2["Name"])]
enough_data = detect_data[["pl_name", "pl_bmassj", "pl_radj", "pl_orbsmax", "sy_dist", "st_age"]]
enough_data = detect_data[["pl_name", "pl_bmassj", "pl_radj", "pl_orbsmax", "sy_dist", "st_age"]]
enough_data.rename(columns={"pl_name": "Name"}, inplace=True)

merged = pd.merge(enough_data, df2, on="Name")


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

def one_format(x):
    return '{:.1f}'.format(x)


def zero_format(x):
    return '{:.0f}'.format(x)

merged['RA'], merged['DEC'] = zip(*(merged["Name"].str[:-1]).apply(query_simbad))
merged["vmax"] = merged["Max. Frequency [MHz]"]
merged.columns = ["Name", "Mass (MJ)", "Radius (RJ)", "a (AU)", "d (pc)", "t (Gyr)", "nupeak (MHz)", "Phi (mJy)", "RA (J2000)", "DEC (J2000)", "vmax"]
final = merged[["Name", "RA (J2000)", "DEC (J2000)", "Mass (MJ)", "Radius (RJ)", "a (AU)", "d (pc)", "t (Gyr)", "vmax", "Phi (mJy)"]]
final = final.sort_values(by="Phi (mJy)", ascending=False)

final["RA (J2000)"] = final["RA (J2000)"].apply(format_ra)
final["DEC (J2000)"] = final["DEC (J2000)"].apply(format_ra)

final["Mass (MJ)"] = final["Mass (MJ)"].apply(three_format)
final["a (AU)"] = (final["a (AU)"]*1e2).apply(one_format)

final["Radius (RJ)"] = final["Radius (RJ)"].apply(two_format)
final["d (pc)"] = final["d (pc)"].apply(zero_format)
final["t (Gyr)"] = final["t (Gyr)"].apply(two_format)
final["vmax"] = final["vmax"].apply(two_format)
final["Phi (mJy)"] = (final["Phi (mJy)"]*1e3).apply(two_format)

final.to_csv("append_table.csv", columns=["Name", "RA (J2000)", "DEC (J2000)", "vmax", "Phi (mJy)"], index=False)
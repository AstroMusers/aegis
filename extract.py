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


def query_simbad(object_name):
    result_table = Simbad.query_object(object_name)
    if result_table is not None:
        ra = result_table['RA'][0]
        dec = result_table['DEC'][0]
        return ra, dec
    else:
        return None, None


merged['RA'], merged['DEC'] = zip(*merged["Name"].apply(query_simbad))


# enough_data["freq"], detect_data["flux"] = freq, flux
# sorted = enough_data.sort_values(by="pl_name")
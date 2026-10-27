"""Shared coordinate and radio-catalog table helpers for Aegis."""

import pandas as pd
from tdpy.astro import format_equatorial_coordinates


def format_coords(row: pd.Series) -> pd.Series:
    """Format decimal-degree coordinates as J2000-style table strings."""

    ra_string, dec_string = format_equatorial_coordinates(row["ra"], row["dec"], precision=0)
    return pd.Series({"RA": ra_string, "DEC": dec_string})


def build_radio_catalog_table(
    exoplanet_catalog: pd.DataFrame,
    predictions: pd.DataFrame,
    frequency_column: str,
    flux_column: str,
) -> pd.DataFrame:
    """Merge radio predictions into the shared, sorted Aegis table schema."""
    planet_columns = [
        "pl_name", "pl_bmassj", "pl_radj", "pl_orbsmax", "sy_dist", "st_age", "ra", "dec"
    ]
    planet_rows = exoplanet_catalog[
        exoplanet_catalog["pl_name"].isin(predictions["Name"])
    ][planet_columns].copy()
    planet_rows.rename(columns={"pl_name": "Name"}, inplace=True)
    prediction_rows = predictions[
        ["Name", frequency_column, flux_column]
    ].rename(columns={frequency_column: "nupeak (MHz)", flux_column: "Phi (mJy)"})

    merged = pd.merge(planet_rows, prediction_rows, on="Name")
    merged[["ra", "dec"]] = merged.apply(format_coords, axis=1)
    merged["vmax"] = merged["nupeak (MHz)"]
    merged.columns = [
        "Name", "Mass (MJ)", "Radius (RJ)", "a (AU)", "d (pc)", "t (Gyr)",
        "RA (J2000)", "DEC (J2000)", "nupeak (MHz)", "Phi (mJy)", "vmax",
    ]
    output_columns = [
        "Name", "RA (J2000)", "DEC (J2000)", "Mass (MJ)", "Radius (RJ)",
        "a (AU)", "d (pc)", "t (Gyr)", "vmax", "Phi (mJy)",
    ]
    return merged[output_columns].sort_values(by="Phi (mJy)", ascending=False)


def convert_numeric_column(column: pd.Series) -> pd.Series:
    """Convert a fully numeric column and retain its original values otherwise."""

    try:
        return pd.to_numeric(column)
    except (TypeError, ValueError):
        return column


def zero_format(x) -> str:
    """Format a numeric value with zero decimal places."""
    return "{:.0f}".format(x)


def one_format(x) -> str:
    """Format a numeric value with one decimal place."""
    return "{:.1f}".format(x)


def two_format(x) -> str:
    """Format a numeric value with two decimal places."""
    return "{:.2f}".format(x)


def three_format(x) -> str:
    """Format a numeric value with three decimal places."""
    return "{:.3f}".format(x)
"""Shared coordinate formatting for Aegis table-generation scripts."""

import pandas as pd

try:
    from tdpy.astro import format_equatorial_coordinates
except ImportError:
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    def format_equatorial_coordinates(ra_deg, dec_deg, precision=0):
        coord = SkyCoord(ra=float(ra_deg) * u.deg, dec=float(dec_deg) * u.deg)
        ra_str = coord.ra.to_string(unit=u.hour, sep=":", precision=precision, pad=True)
        dec_str = coord.dec.to_string(unit=u.deg, sep=":", precision=precision, alwayssign=True, pad=True)
        return ra_str, dec_str


def format_coords(row: pd.Series) -> pd.Series:
    """Format decimal-degree coordinates as J2000-style table strings."""

    ra_string, dec_string = format_equatorial_coordinates(row["ra"], row["dec"], precision=0)
    return pd.Series({"RA": ra_string, "DEC": dec_string})


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
import numpy as np
import pandas as pd
import pytest

from coordinate_formatting import (
    build_radio_catalog_table,
    convert_numeric_column,
    format_coords,
    one_format,
    three_format,
    two_format,
    zero_format,
)


def test_format_coords_preserves_table_format():
    """Decimal-degree coordinates retain the existing table representation."""

    formatted = format_coords(pd.Series({"ra": 0.0, "dec": -30.0}))

    assert formatted.to_dict() == {"RA": "00:00:00", "DEC": "-30:00:00"}


def test_coordinate_formatting_uses_tdpy_owner():
    from coordinate_formatting import format_equatorial_coordinates
    from tdpy.astro import format_equatorial_coordinates as tdpy_formatter

    assert format_equatorial_coordinates is tdpy_formatter


@pytest.mark.parametrize(
    ("frequency_column", "flux_column"),
    [("Freq(MHz)", "Flux(mJy)"), ("Max. Frequency [MHz]", "Max. Flux Density [Jy]")],
)
def test_build_radio_catalog_table_preserves_schema_sort_and_inputs(
    frequency_column, flux_column
):
    planets = pd.DataFrame([
        {
            "pl_name": "Planet A",
            "pl_bmassj": 1.0,
            "pl_radj": 1.0,
            "pl_orbsmax": 0.1,
            "sy_dist": 10.0,
            "st_age": 5.0,
            "ra": 180.0,
            "dec": 20.0,
        },
        {
            "pl_name": "Planet B",
            "pl_bmassj": 2.0,
            "pl_radj": 1.5,
            "pl_orbsmax": 0.2,
            "sy_dist": 20.0,
            "st_age": 6.0,
            "ra": 0.0,
            "dec": -30.0,
        },
    ])
    predictions = pd.DataFrame([
        {"Name": "Planet B", frequency_column: 10.0, flux_column: 2.0},
        {"Name": "Planet A", frequency_column: 40.0, flux_column: 5.0},
    ])
    original_planets = planets.copy(deep=True)
    original_predictions = predictions.copy(deep=True)

    table = build_radio_catalog_table(
        planets, predictions, frequency_column, flux_column
    )

    assert table["Name"].tolist() == ["Planet A", "Planet B"]
    assert table["RA (J2000)"].tolist() == ["12:00:00", "00:00:00"]
    assert table["DEC (J2000)"].tolist() == ["+20:00:00", "-30:00:00"]
    assert table["vmax"].tolist() == [40.0, 10.0]
    assert table["Phi (mJy)"].tolist() == [5.0, 2.0]
    pd.testing.assert_frame_equal(planets, original_planets)
    pd.testing.assert_frame_equal(predictions, original_predictions)


def test_numeric_precision_formatters():
    assert zero_format(3.14159) == "3"
    assert one_format(3.14159) == "3.1"
    assert two_format(3.14159) == "3.14"
    assert three_format(3.14159) == "3.142"


def test_convert_numeric_column_preserves_nonnumeric_values():
    numeric = convert_numeric_column(pd.Series(["1", "2"]))
    mixed = convert_numeric_column(pd.Series(["1", "text"]))

    assert numeric.tolist() == [1, 2]
    assert mixed.tolist() == ["1", "text"]


def test_default_calculate_returns_positive_emission():
    from does_it_make_a_difference import default_calculate

    freq_mhz, flux_jy = default_calculate("tau Boo b")
    assert freq_mhz > 0.0
    assert flux_jy > 0.0


def test_keplerian_speed_uses_tdpy_implementation_for_arrays():
    from radio_module import keplerian
    from tdpy.exoplanet import keplerian_orbital_speed

    assert keplerian is keplerian_orbital_speed
    np.testing.assert_array_equal(
        keplerian(np.array([1.0, 2.0]), np.array([1.0, 4.0])),
        keplerian_orbital_speed(np.array([1.0, 2.0]), np.array([1.0, 4.0])),
    )


def test_sketch_and_spiral_plots_execute_without_display():
    from sketch import draw_magnetosphere_sketch
    from parker_spiral import plot_parker_spiral
    from radius_valley import plot_radius_valleys
    from visibility import max_elev, time_above_elevation, visibility_plot

    draw_magnetosphere_sketch(output_path=None, show=False)
    plot_parker_spiral(save=False, show=False)
    plot_radius_valleys(save=False, show=False)

    elev = max_elev(0.0, 30.0, 54.9)
    assert elev > 0.0

    hours = time_above_elevation(0.0, 30.0, 54.9, 20.0)
    assert hours > 0.0

    visibility_plot(save=False, show=False)

import pandas as pd

from coordinate_formatting import (
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


def test_numeric_precision_formatters():
    assert zero_format(3.14159) == "3"
    assert one_format(3.14159) == "3.1"
    assert two_format(3.14159) == "3.14"
    assert three_format(3.14159) == "3.142"


def test_default_calculate_returns_positive_emission():
    from does_it_make_a_difference import default_calculate

    freq_mhz, flux_jy = default_calculate("tau Boo b")
    assert freq_mhz > 0.0
    assert flux_jy > 0.0


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

import pandas as pd
import pytest
from tdpy.verbosity import print

import append


@pytest.mark.parametrize("verbosity", ["0", "1"])
def test_append_prediction_file_narration_respects_verbosity(
    tmp_path, monkeypatch, capsys, verbosity
):
    monkeypatch.setenv("TDPY_VERBOSITY", verbosity)
    monkeypatch.setattr(append, "DATA_DIR", tmp_path)
    output_tables = tmp_path / "Output Tables"
    output_tables.mkdir()
    nasa_path = tmp_path / "NASA3004.csv"
    predictions_path = output_tables / "all.csv"
    print(f"Writing to {nasa_path}...")
    pd.DataFrame(
        [
            {
                "pl_name": "Planet b",
                "pl_bmassj": 1.0,
                "pl_radj": 1.0,
                "pl_orbsmax": 0.1,
                "sy_dist": 10.0,
                "st_age": 5.0,
                "ra": 180.0,
                "dec": 20.0,
            }
        ]
    ).to_csv(nasa_path, index=False)
    print(f"Writing to {predictions_path}...")
    pd.DataFrame(
        [
            {"Name": "Planet b", "Max. Frequency [MHz]": 2.5, "Max. Flux Density [Jy]": 4.0},
            {"Name": "Other b", "Max. Frequency [MHz]": 11.0, "Max. Flux Density [Jy]": 1.0},
        ]
    ).to_csv(predictions_path, index=False)
    capsys.readouterr()

    final = append.append_predictions(save=True)

    assert len(final) == 1
    filtered_path = output_tables / "0.1-10MHz.csv"
    output_path = tmp_path / "append_table.csv"
    assert filtered_path.is_file()
    assert output_path.is_file()
    if verbosity == "0":
        assert capsys.readouterr().out == ""
    else:
        assert capsys.readouterr().out.splitlines() == [
            f"Reading from {nasa_path}...",
            f"Reading from {predictions_path}...",
            "2",
            f"Writing to {filtered_path}...",
            f"Writing to {output_path}...",
        ]
import csv
from pathlib import Path

import pytest
from tdpy.verbosity import print

import extract


@pytest.mark.parametrize("verbosity", ["0", "1"])
def test_extract_file_narration_respects_verbosity(tmp_path, monkeypatch, capsys, verbosity):
    monkeypatch.setenv("TDPY_VERBOSITY", verbosity)
    monkeypatch.setattr(extract, "DATA_DIR", tmp_path)
    input_rows = {
        "NASA3004.csv": [
            {
                "pl_name": "Planet b",
                "pl_bmassj": "1.0",
                "pl_radj": "1.0",
                "pl_orbsmax": "0.1",
                "sy_dist": "10.0",
                "st_age": "5.0",
                "ra": "180.0",
                "dec": "20.0",
            }
        ],
        "detectables_both.csv": [
            {"Name": "Planet b", "Freq(MHz)": "40.0", "Flux(mJy)": "2.0"}
        ],
    }
    for filename, rows in input_rows.items():
        path = tmp_path / filename
        print(f"Writing to {path}...")
        with path.open("w", newline="", encoding="utf-8") as output_file:
            writer = csv.DictWriter(output_file, fieldnames=rows[0])
            writer.writeheader()
            writer.writerows(rows)
    capsys.readouterr()

    final, more_info = extract.extract_observable_targets(save=True)

    assert len(final) == len(more_info) == 1
    output_path = tmp_path / "obs_table.csv"
    assert output_path.is_file()
    if verbosity == "0":
        assert capsys.readouterr().out == ""
    else:
        assert capsys.readouterr().out.splitlines() == [
            f"Reading from {tmp_path / 'NASA3004.csv'}...",
            f"Reading from {tmp_path / 'detectables_both.csv'}...",
            f"Writing to {output_path}...",
        ]
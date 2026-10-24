import csv
import json
from pathlib import Path

import pytest

from publish_radio_catalog import CATALOG_COLUMNS, publish_radio_catalog
from tdpy.verbosity import print


def write_catalog(path: Path, rows: list[dict[str, str]]) -> None:
    print(f"Writing to {path}...")
    with path.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=CATALOG_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def test_publisher_writes_csv_and_json_for_catalog(tmp_path: Path) -> None:
    source_path = tmp_path / "append_table.csv"
    output_directory = tmp_path / "docs"
    rows = [
        {
            "Name": "Planet b",
            "RA (J2000)": "12:00:00",
            "DEC (J2000)": "+20:00:00",
            "vmax": "2.5",
            "Phi (Jy)": "4.0",
        }
    ]
    write_catalog(source_path, rows)

    csv_path, json_path = publish_radio_catalog(
        source_path,
        output_directory,
        boolcheckphys=True,
    )

    print(f"Reading from {csv_path}...")
    with csv_path.open(newline="", encoding="utf-8") as input_file:
        written_rows = list(csv.DictReader(input_file))
    print(f"Reading from {json_path}...")
    with json_path.open(encoding="utf-8") as input_file:
        json_rows = json.load(input_file)
    assert written_rows == rows
    assert json_rows == rows


@pytest.mark.parametrize("verbosity, expected_output", [("0", ""), ("1", None)])
def test_publisher_file_narration_respects_verbosity(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    verbosity: str,
    expected_output: str | None,
) -> None:
    source_path = tmp_path / "append_table.csv"
    output_directory = tmp_path / "docs"
    monkeypatch.setenv("TDPY_VERBOSITY", verbosity)
    write_catalog(
        source_path,
        [{
            "Name": "Planet b",
            "RA (J2000)": "12:00:00",
            "DEC (J2000)": "+20:00:00",
            "vmax": "2.5",
            "Phi (Jy)": "4.0",
        }],
    )
    capsys.readouterr()

    csv_path, json_path = publish_radio_catalog(source_path, output_directory)

    if expected_output is None:
        expected_output = (
            f"Reading from {source_path}...\n"
            f"Writing to {csv_path}...\n"
            f"Writing to {json_path}...\n"
        )
    assert capsys.readouterr().out == expected_output


def test_publisher_rejects_missing_catalog_columns(tmp_path: Path) -> None:
    source_path = tmp_path / "incomplete.csv"
    print(f"Writing to {source_path}...")
    source_path.write_text("Name,vmax\nPlanet b,2.5\n", encoding="utf-8")

    with pytest.raises(ValueError, match="missing required columns"):
        publish_radio_catalog(source_path, tmp_path / "docs")


@pytest.mark.parametrize(
    "frequency, flux_density",
    [("0", "4"), ("nan", "4"), ("2.5", "0"), ("2.5", "inf")],
)
def test_publisher_optional_physicality_check_rejects_invalid_predictions(
    tmp_path: Path,
    frequency: str,
    flux_density: str,
) -> None:
    source_path = tmp_path / "invalid.csv"
    write_catalog(
        source_path,
        [{
            "Name": "Planet b",
            "RA (J2000)": "12:00:00",
            "DEC (J2000)": "+20:00:00",
            "vmax": frequency,
            "Phi (Jy)": flux_density,
        }],
    )

    with pytest.raises(ValueError):
        publish_radio_catalog(source_path, tmp_path / "docs", boolcheckphys=True)


def test_catalog_page_exposes_accessible_filters_and_sorting() -> None:
    page_path = Path(__file__).parent / "docs" / "index.html"
    print(f"Reading from {page_path}...")
    page = page_path.read_text(encoding="utf-8")

    assert 'id="search"' in page
    assert 'aria-live="polite"' in page
    assert 'data-sort="vmax"' in page
    assert 'data-sort="Phi (Jy)"' in page
    assert 'fetch("radio_catalog.json")' in page
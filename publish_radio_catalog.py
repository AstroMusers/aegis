"""Publish the Aegis predicted-radio-target table as static site data."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

from tdpy.verbosity import print

from repository_paths import get_repository_path


CATALOG_COLUMNS = ("Name", "RA (J2000)", "DEC (J2000)", "vmax", "Phi (Jy)")


def publish_radio_catalog(
    source_path: str | Path | None = None,
    output_directory: str | Path | None = None,
    *,
    boolcheckphys: bool = False,
) -> tuple[Path, Path]:
    """Write catalog CSV and JSON files for the static website."""
    repository_path = get_repository_path()
    source_path = Path(source_path) if source_path is not None else repository_path / "append_table.csv"
    output_directory = (
        Path(output_directory) if output_directory is not None else repository_path / "docs"
    )

    print(f"Reading from {source_path}...")
    with source_path.open(newline="", encoding="utf-8") as source_file:
        reader = csv.DictReader(source_file)
        missing_columns = set(CATALOG_COLUMNS).difference(reader.fieldnames or ())
        if missing_columns:
            raise ValueError(f"Catalog is missing required columns: {sorted(missing_columns)}")
        rows = [{column: row[column] for column in CATALOG_COLUMNS} for row in reader]

    if not rows:
        raise ValueError("Catalog must contain at least one target")
    if boolcheckphys:
        for row in rows:
            frequency_mhz = float(row["vmax"])
            flux_density_jy = float(row["Phi (Jy)"])
            if not math.isfinite(frequency_mhz) or frequency_mhz <= 0.0:
                raise ValueError(f"Target {row['Name']!r} has an invalid maximum frequency")
            if not math.isfinite(flux_density_jy) or flux_density_jy <= 0.0:
                raise ValueError(f"Target {row['Name']!r} has an invalid predicted flux density")

    output_directory.mkdir(parents=True, exist_ok=True)
    csv_path = output_directory / "radio_catalog.csv"
    json_path = output_directory / "radio_catalog.json"

    print(f"Writing to {csv_path}...")
    with csv_path.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=CATALOG_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Writing to {json_path}...")
    with json_path.open("w", encoding="utf-8") as output_file:
        json.dump(rows, output_file, ensure_ascii=False, indent=2)
        output_file.write("\n")
    return csv_path, json_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check-physicality",
        action="store_true",
        help="check predicted frequency and flux values before publishing",
    )
    arguments = parser.parse_args()
    publish_radio_catalog(boolcheckphys=arguments.check_physicality)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
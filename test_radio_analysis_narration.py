import os
import subprocess
import sys
from pathlib import Path

import pytest


@pytest.mark.parametrize("verbosity", ["0", "1"])
def test_radio_analysis_import_narrates_input_files(verbosity):
    repository_path = Path(__file__).resolve().parent
    environment = os.environ.copy()
    python_paths = [str(repository_path)]
    python_paths.extend(
        path for path in environment.get("PYTHONPATH", "").split(os.pathsep) if path
    )
    environment.update(
        {
            "AEGIS_PATH": str(repository_path),
            "PYTHONPATH": os.pathsep.join(python_paths),
            "TDPY_VERBOSITY": verbosity,
        }
    )
    result = subprocess.run(
        [sys.executable, "-c", "import does_it_make_a_difference"],
        cwd=repository_path,
        env=environment,
        capture_output=True,
        check=True,
        text=True,
    )

    if verbosity == "0":
        assert result.stdout == ""
    else:
        assert result.stdout.splitlines() == [
            f"Reading from {repository_path / 'rotation.csv'}...",
            f"Reading from {repository_path / 'rotation.csv'}...",
            f"Reading from {repository_path / 'NASA3004.csv'}...",
            f"Reading from {repository_path / 'stellar_wind_predictions_3004.tsv'}...",
        ]
import os
import subprocess
import sys
from pathlib import Path

import pytest


@pytest.mark.parametrize("verbosity", ["0", "1"])
def test_visibility_file_narration_respects_verbosity(tmp_path, verbosity):
    repository_path = Path(__file__).resolve().parent
    output_path = tmp_path / f"final-table-{verbosity}.csv"
    environment = os.environ.copy()
    environment.update(
        {
            "AEGIS_PATH": str(repository_path),
            "PYTHONPATH": str(repository_path),
            "TDPY_VERBOSITY": verbosity,
            "VISIBILITY_OUTPUT": str(output_path),
        }
    )
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "import os; import visibility; "
            "visibility.generate_final_table(output_path=os.environ['VISIBILITY_OUTPUT'])",
        ],
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
            f"Reading from {repository_path / 'obs_table.csv'}...",
            f"Reading from {repository_path / 'observables.npz'}...",
            f"Reading from {repository_path / 'obs_table.csv'}...",
            f"Writing to {output_path}...",
        ]
    assert output_path.is_file()
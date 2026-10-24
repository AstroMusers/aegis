import os
import subprocess
import sys
from pathlib import Path

import pytest

@pytest.mark.parametrize("verbosity", ["0", "1"])
def test_plotter_import_narrates_input_files(verbosity):
    repository_path = Path(__file__).resolve().parent
    environment = os.environ.copy()
    python_paths = [str(repository_path)]
    if environment.get("PYTHONPATH"):
        python_paths.append(environment["PYTHONPATH"])
    environment.update(
        {
            "AEGIS_PATH": str(repository_path),
            "PYTHONPATH": os.pathsep.join(python_paths),
            "TDPY_VERBOSITY": verbosity,
        }
    )
    import_code = (
        "import pandas as pd; import plotter; from coordinate_formatting import convert_numeric_column; "
        "assert convert_numeric_column(pd.Series(['1', '2'])).tolist() == [1, 2]; "
        "assert convert_numeric_column(pd.Series(['x', '2'])).tolist() == ['x', '2']"
    )
    result = subprocess.run(
        [sys.executable, "-c", import_code],
        cwd=repository_path,
        env=environment,
        capture_output=True,
        check=False,
        text=True,
    )
    assert result.returncode == 0, result.stderr

    if verbosity == "0":
        assert result.stdout == ""
    else:
        assert result.stdout.splitlines() == [
            f"Reading from {repository_path / 'sensitivities.csv'}...",
            f"Reading from {repository_path / 'nenufar.npz'}...",
            f"Reading from {repository_path / 'df1.csv'}...",
            f"Reading from {repository_path / 'all.npz'}...",
        ]
import os
import subprocess
import sys
from pathlib import Path

import pytest


@pytest.mark.parametrize("verbosity", ["0", "1"])
def test_radius_valley_file_narration_respects_verbosity(verbosity):
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
    code = (
        "import matplotlib.pyplot as plt; import radius_valley; "
        "plt.savefig = lambda *args, **kwargs: None; "
        "radius_valley.plot_radius_valleys(save=True)"
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
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
            f"Reading from {repository_path / 'NASA3004.csv'}...",
            f"Reading from {repository_path / 'Output Tables' / 'all.csv'}...",
            f"Writing to {repository_path / 'valleys.pdf'}...",
            f"Writing to {repository_path / 'single_valley.pdf'}...",
        ]
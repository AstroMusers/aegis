import importlib
from pathlib import Path

import numpy as np
import pytest
from tdpy.verbosity import print

import wind


@pytest.mark.parametrize("verbosity", ["0", "1"])
def test_wind_pipeline_file_narration_respects_verbosity(
    tmp_path, monkeypatch, capsys, verbosity
):
    repository_path = Path(__file__).resolve().parent
    monkeypatch.setenv("AEGIS_PATH", str(repository_path))
    monkeypatch.setenv("TDPY_VERBOSITY", verbosity)
    module = importlib.reload(wind)
    import_output = capsys.readouterr().out

    monkeypatch.setattr(module, "DATA_DIR", tmp_path)
    monkeypatch.setattr(module, "filename", str(tmp_path / "NASA3004.csv"))
    monkeypatch.setattr(module, "names", np.array(["tau Boo b"]))
    monkeypatch.setattr(module, "semis", np.array([0.1]))
    monkeypatch.setattr(module, "masses", np.array([1.0]))
    monkeypatch.setattr(module, "ages", np.array([5.0]))
    monkeypatch.setattr(module, "specs", np.array(["G"]))
    monkeypatch.setattr(module, "v_at_1_au", lambda age: 1.0)
    monkeypatch.setattr(module, "critical_radius", lambda mass, temperature: 1.0)
    monkeypatch.setattr(module, "fsolve", lambda function, guess: np.asarray(guess))
    monkeypatch.setattr(module, "non_decreasing", lambda profile: True)

    module.run_wind_simulation(save=True, show=False)
    run_output = capsys.readouterr().out

    assert (tmp_path / "winds.pdf").is_file()
    assert (tmp_path / "stellar_wind_predictions_3004.tsv").is_file()
    assert (tmp_path / "taub_wind.npz").is_file()
    if verbosity == "0":
        assert import_output == run_output == ""
    else:
        assert import_output.splitlines() == [
            f"Reading from {repository_path / 'NASA3004.csv'}..."
        ]
        assert run_output.splitlines() == [
            f"Writing to {tmp_path / 'winds.pdf'}...",
            "All good, fantastic work!",
            f"Writing to {tmp_path / 'stellar_wind_predictions_3004.tsv'}...",
            f"Writing to {tmp_path / 'taub_wind.npz'}...",
        ]
import numpy as np
import pytest
from tdpy.verbosity import print

import parker_spiral


@pytest.mark.parametrize("verbosity", ["0", "1"])
def test_parker_spiral_file_narration_respects_verbosity(
    tmp_path, monkeypatch, capsys, verbosity
):
    monkeypatch.setenv("TDPY_VERBOSITY", verbosity)
    monkeypatch.setattr(parker_spiral, "DATA_DIR", tmp_path)
    monkeypatch.setattr(
        parker_spiral.np,
        "load",
        lambda path: {
            "ranges": np.logspace(-3, 0, 50),
            "speeds": np.full(50, 100.0),
        },
    )

    parker_spiral.plot_parker_spiral(save=True)

    input_path = tmp_path / "taub_wind.npz"
    output_path = tmp_path / "parker_spiral.pdf"
    assert output_path.is_file()
    assert output_path.read_bytes().startswith(b"%PDF-")
    if verbosity == "0":
        assert capsys.readouterr().out == ""
    else:
        assert capsys.readouterr().out.splitlines() == [
            f"Reading from {input_path}...",
            f"Writing to {output_path}...",
        ]
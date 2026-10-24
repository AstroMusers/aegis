import pytest
from tdpy.verbosity import print

from sketch import draw_magnetosphere_sketch


@pytest.mark.parametrize("verbosity", ["0", "1"])
def test_sketch_file_narration_respects_verbosity(tmp_path, monkeypatch, capsys, verbosity):
    monkeypatch.setenv("TDPY_VERBOSITY", verbosity)
    output_path = tmp_path / "magnetosphere.png"

    draw_magnetosphere_sketch(output_path=output_path)

    assert output_path.is_file()
    print(f"Reading from {output_path}...")
    assert output_path.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    if verbosity == "0":
        assert capsys.readouterr().out == ""
    else:
        assert capsys.readouterr().out.splitlines() == [
            f"Writing to {output_path}...",
            f"Reading from {output_path}...",
        ]
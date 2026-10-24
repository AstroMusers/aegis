import importlib
import os
from pathlib import Path

from packaging.requirements import Requirement

import repository_paths


def test_repository_path_defaults_to_aegis_checkout(monkeypatch):
    monkeypatch.delenv("AEGIS_PATH", raising=False)
    module = importlib.reload(repository_paths)

    expected_path = Path(module.__file__).resolve().parent
    assert os.environ["AEGIS_PATH"] == str(expected_path)
    assert module.get_repository_path() == expected_path
    assert module.get_data_path() == expected_path / "data"
    assert module.get_visuals_path() == expected_path / "visuals"


def test_repository_path_honors_configured_override(tmp_path, monkeypatch):
    monkeypatch.setenv("AEGIS_PATH", str(tmp_path))

    assert repository_paths.get_repository_path() == tmp_path.resolve()


def test_requirements_are_valid_and_use_no_local_file_urls():
    requirements_path = Path(__file__).resolve().parent / "requirements.txt"
    requirements = [
        Requirement(line)
        for line in requirements_path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]

    assert requirements
    assert all(
        requirement.url is None or not requirement.url.startswith("file:")
        for requirement in requirements
    )

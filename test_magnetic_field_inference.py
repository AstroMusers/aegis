import json
from pathlib import Path
import sys
from types import ModuleType

import numpy as np
import pytest

from magnetic_field_inference import _cutoff_log_likelihood, infer_surface_field


def test_cutoff_notebook_uses_pcat_inference_and_keeps_paths_private():
    path = Path(__file__).parent / "examples" / "magnetic_field_cutoff.ipynb"
    print("Reading from examples/magnetic_field_cutoff.ipynb...")
    notebook = json.loads(path.read_text(encoding="utf-8"))
    code = "\n".join(
        "".join(cell["source"])
        for cell in notebook["cells"]
        if cell["cell_type"] == "code"
    )
    document = json.dumps(notebook)

    assert "from magnetic_field_inference import infer_surface_field" in code
    assert "infer_surface_field(" in code
    assert "PCAT" in document
    assert "not telescope measurements" in document
    assert "/" + "Users/" not in document
    assert 'repository_root / "examples" / "visuals"' in code
    figure_path = Path(__file__).parent / "examples" / "visuals" / "aegis_synthetic_cutoff_field_posterior.png"
    assert figure_path.is_file() and figure_path.stat().st_size > 0


def test_cutoff_likelihood_peaks_at_aegis_field_relation():
    observed_frequency_mhz = 40.0 * 24.0 / 14.0
    state = {"observed_frequency_mhz": observed_frequency_mhz, "uncertainty_mhz": 1.0}

    at_truth = _cutoff_log_likelihood([24.0], state)
    away_from_truth = _cutoff_log_likelihood([48.0], state)

    assert at_truth > away_from_truth
    assert np.isfinite(at_truth)


@pytest.mark.parametrize("field_gauss", [12.0, 24.0, 48.0])
@pytest.mark.parametrize("uncertainty_mhz", [0.02, 1.0, 9.0])
def test_cutoff_likelihood_preserves_normalized_measurement_density(field_gauss, uncertainty_mhz):
    from radio_module import max_freq

    observed = 40.0 * 24.0 / 14.0  # [MHz]
    predicted = max_freq(field_gauss) / 1.0e6  # [MHz]
    residual = (observed - predicted) / uncertainty_mhz
    expected = -0.5 * (residual**2 + np.log(2.0 * np.pi * uncertainty_mhz**2))
    state = {"observed_frequency_mhz": observed, "uncertainty_mhz": uncertainty_mhz}
    assert _cutoff_log_likelihood([field_gauss], state) == pytest.approx(expected, rel=1e-14)


@pytest.mark.parametrize("field_gauss", [0.0, -1.0, np.nan, np.inf])
def test_cutoff_likelihood_retains_invalid_field_rejection(field_gauss):
    assert _cutoff_log_likelihood([field_gauss], {}) == -np.inf


def test_infer_surface_field_delegates_to_pcat(monkeypatch):
    captured = {}
    fake_pcat = ModuleType("pcat")
    fake_pcat.__path__ = []
    fake_fixed = ModuleType("pcat.fixed")

    def sample_fixed_chains(*args, **kwargs):
        captured["args"] = args
        captured["kwargs"] = kwargs
        return np.full((2, 5, 1), 24.0), np.zeros((2, 5))

    fake_fixed.sample_fixed_chains = sample_fixed_chains
    monkeypatch.setitem(sys.modules, "pcat", fake_pcat)
    monkeypatch.setitem(sys.modules, "pcat.fixed", fake_fixed)

    posterior = infer_surface_field(
        40.0 * 24.0 / 14.0,
        1.0,
        field_bounds_gauss=(1.0, 100.0),
        chain_count=2,
        sample_count=10,
        burn_count=5,
    )

    assert captured["args"][3] == ("surface_field_gauss",)
    assert captured["kwargs"]["typeseed"] == 0
    assert posterior.chain_gauss.shape == (2, 5)
    assert posterior.field_interval_gauss == pytest.approx((24.0, 24.0, 24.0))


def test_infer_surface_field_rejects_invalid_measurement():
    with pytest.raises(ValueError, match="uncertainty_mhz"):
        infer_surface_field(50.0, 0.0)
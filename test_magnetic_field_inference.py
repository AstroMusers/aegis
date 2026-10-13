import sys
from types import ModuleType

import numpy as np
import pytest

from magnetic_field_inference import _cutoff_log_likelihood, infer_surface_field


def test_cutoff_likelihood_peaks_at_aegis_field_relation():
    observed_frequency_mhz = 40.0 * 24.0 / 14.0
    state = {"observed_frequency_mhz": observed_frequency_mhz, "uncertainty_mhz": 1.0}

    at_truth = _cutoff_log_likelihood([24.0], state)
    away_from_truth = _cutoff_log_likelihood([48.0], state)

    assert at_truth > away_from_truth
    assert np.isfinite(at_truth)


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
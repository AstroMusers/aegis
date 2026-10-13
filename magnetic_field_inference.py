"""Magnetic-field inference from Aegis radio-cutoff measurements."""

from dataclasses import dataclass

import numpy as np

from radio_module import max_freq


@dataclass(frozen=True)
class MagneticFieldPosterior:
    """Posterior samples for a planetary surface field inferred from a cutoff."""

    chain_gauss: np.ndarray
    log_probability: np.ndarray
    observed_frequency_mhz: float
    uncertainty_mhz: float
    prior_bounds_gauss: tuple[float, float]

    @property
    def field_interval_gauss(self) -> tuple[float, float, float]:
        """Return median and central 68% interval [G]."""

        quantiles = np.percentile(self.chain_gauss, [16.0, 50.0, 84.0])
        return tuple(float(value) for value in quantiles)


def _cutoff_log_likelihood(parameters, state) -> float:
    field_gauss = float(parameters[0])  # [G]
    if not np.isfinite(field_gauss) or field_gauss <= 0.0:
        return -np.inf
    predicted_frequency_mhz = float(max_freq(field_gauss) / 1.0e6)  # [MHz]
    residual = (state["observed_frequency_mhz"] - predicted_frequency_mhz) / state["uncertainty_mhz"]
    return float(-0.5 * (residual**2 + np.log(2.0 * np.pi * state["uncertainty_mhz"] ** 2)))


def infer_surface_field(
    observed_cutoff_frequency_mhz: float,
    uncertainty_mhz: float,
    *,
    field_bounds_gauss: tuple[float, float] = (0.1, 200.0),
    chain_count: int = 4,
    sample_count: int = 1200,
    burn_count: int = 600,
    seed: int = 0,
    pathbase: str | None = None,
) -> MagneticFieldPosterior:
    """Infer planetary surface field from a detected CMI cutoff frequency.

    The fit uses Aegis's cyclotron-cutoff relation, a Gaussian measurement
    likelihood, and a uniform field prior. It complements the forward Monte
    Carlo uncertainty propagation with an observation-driven parameter estimate.
    """

    if not np.isfinite(observed_cutoff_frequency_mhz) or observed_cutoff_frequency_mhz <= 0.0:
        raise ValueError("observed_cutoff_frequency_mhz must be finite and positive")
    if not np.isfinite(uncertainty_mhz) or uncertainty_mhz <= 0.0:
        raise ValueError("uncertainty_mhz must be finite and positive")
    lower_field, upper_field = field_bounds_gauss
    if not np.isfinite([lower_field, upper_field]).all() or lower_field <= 0.0 or lower_field >= upper_field:
        raise ValueError("field_bounds_gauss must be finite, positive, and increasing")
    if chain_count < 1 or sample_count <= burn_count or burn_count < 0:
        raise ValueError("chain_count must be positive and sample_count must exceed burn_count")

    try:
        from pcat.fixed import sample_fixed_chains
    except ImportError as exception:
        raise ImportError(
            "Aegis Bayesian inference dependencies are missing; install requirements.txt."
        ) from exception

    state = {
        "observed_frequency_mhz": float(observed_cutoff_frequency_mhz),
        "uncertainty_mhz": float(uncertainty_mhz),
    }
    chains, log_probability = sample_fixed_chains(
        state,
        _cutoff_log_likelihood,
        None,
        ("surface_field_gauss",),
        ("self",),
        np.asarray([lower_field]),
        np.asarray([upper_field]),
        None,
        None,
        None,
        chain_count,
        sample_count,
        burn_count,
        pathbase=pathbase,
        typeverb=0,
        typeseed=seed,
        boolmakeplot=False,
        boolmakeplotinit=False,
        makeanim=False,
        booldiag=False,
    )
    return MagneticFieldPosterior(
        chain_gauss=np.asarray(chains, dtype=float)[:, :, 0],
        log_probability=np.asarray(log_probability, dtype=float),
        observed_frequency_mhz=float(observed_cutoff_frequency_mhz),
        uncertainty_mhz=float(uncertainty_mhz),
        prior_bounds_gauss=(float(lower_field), float(upper_field)),
    )
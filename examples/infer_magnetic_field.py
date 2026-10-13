#!/usr/bin/env python3
"""Infer a field from an explicitly simulated cyclotron cutoff measurement."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from radio_module import max_freq
from magnetic_field_inference import infer_surface_field


DATA_DIR = Path(__file__).resolve().parents[1]
TRUE_FIELD_GAUSS = 24.0  # [G], assumed synthetic truth
MEASUREMENT_UNCERTAINTY_MHZ = 2.0  # [MHz], assumed synthetic uncertainty


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-count", type=int, default=1200)
    parser.add_argument("--burn-count", type=int, default=600)
    parser.add_argument("--chain-count", type=int, default=4)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output", type=Path, default=DATA_DIR / "visuals" / "magnetic_field_posterior.png")
    arguments = parser.parse_args()

    rng = np.random.default_rng(arguments.seed)
    true_frequency_mhz = max_freq(TRUE_FIELD_GAUSS) / 1.0e6  # [MHz]
    observed_frequency_mhz = float(
        rng.normal(true_frequency_mhz, MEASUREMENT_UNCERTAINTY_MHZ)
    )  # [MHz]
    posterior = infer_surface_field(
        observed_frequency_mhz,
        MEASUREMENT_UNCERTAINTY_MHZ,
        chain_count=arguments.chain_count,
        sample_count=arguments.sample_count,
        burn_count=arguments.burn_count,
        seed=arguments.seed,
        pathbase=str(DATA_DIR / "data" / "pcat_cutoff_field"),
    )
    lower_field, median_field, upper_field = posterior.field_interval_gauss
    print(
        f"Synthetic cutoff: {observed_frequency_mhz:.2f} +/- "
        f"{MEASUREMENT_UNCERTAINTY_MHZ:.2f} MHz; inferred field: "
        f"{median_field:.2f} (+{upper_field - median_field:.2f}/"
        f"-{median_field - lower_field:.2f}) G; assumed truth: {TRUE_FIELD_GAUSS:.1f} G."
    )

    fields = posterior.chain_gauss.ravel()
    figure, axis = plt.subplots(figsize=(6.5, 3.8))
    axis.hist(fields, bins=35, density=True, color="#007360", alpha=0.78)
    axis.axvline(TRUE_FIELD_GAUSS, color="#A51417", ls="--", lw=1.2,
                 label="Assumed synthetic truth")
    axis.axvline(median_field, color="#222222", lw=1.1,
                 label="Posterior median")
    axis.set(xlabel="Planetary surface magnetic field [G]",
             ylabel="Posterior density [G$^{-1}$]")
    axis.legend(loc="upper right", fontsize=8)
    axis.grid(False)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    print(f"Writing to {arguments.output}...")
    figure.savefig(arguments.output, dpi=250, bbox_inches="tight", facecolor="white")
    plt.close(figure)


if __name__ == "__main__":
    main()
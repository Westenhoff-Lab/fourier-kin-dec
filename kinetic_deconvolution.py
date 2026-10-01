#!/usr/bin/env python3
"""
Kinetic mode decomposition of time-resolved difference structure factors.

Given:
  - a set of per-timepoint difference-amplitude files (.phs: h k l F weight phase)
  - matching per-timepoint sigma files (.hkl: h k l F sigma)
  - a dark-state sigma file (.hkl)
  - a concentration matrix (CSV) giving the fractional occupancy of each
    kinetic species at each timepoint

this script solves, by linear least squares,

    C @ M = D

for the pure "kinetic mode" difference structure factors M of each species,
where C is the concentration matrix (n_timepoints x n_species) and D is the
stacked difference-amplitude matrix (n_timepoints x n_reflections).

Sigmas are propagated per species using the concentration weights, and used
to build a per-reflection weight for each output kinetic-mode .phs file.

Two files are written per species:
  - a .phs file (h,k,l,Amplitude,weight,Phase) with difference amplitudes, weights and phases for calculating maps and refinement
  - a .hkl file (h,k,l,Amplitude,Sigma) with the difference amplitudes and propagated sigmas

For experimental (as opposed to simulated) data, numerical instability in the
least-squares solve can occasionally produce NaNs or unphysically large
amplitudes/sigmas for a small number of reflections. Any voxel with a NaN, or with
an amplitude or sigma above --amplitude-cutoff, is dropped from both output
files for that species.

Example
-------
    python kinetic_deconvolution.py \\
        --concentrations concentrations.csv \\
        --phs-dir data/phs --phs-pattern "{label}.phs" \\
        --hkl-dir data/hkl --hkl-pattern "{label}_sigma.hkl" \\
        --dark-sigma dark_scaled.hkl \\
        --dark-phase dark_phase.hkl \\
        --min-occurrence 12 \\
        --out-prefix kinetic_mode

concentrations.csv layout (one row per timepoint, in processing order):

    timepoint,pR0,pR1,pR2,pB0
    100ps,0.080,0.010,0.000,0.000
    300ps,0.060,0.025,0.000,0.000
    ...

The first column is the timepoint label used to locate the matching .phs/.hkl
files (via --phs-pattern / --hkl-pattern); the remaining columns are the
per-species occupancies and also become the output-file species names.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from collections import Counter

PHS_COLUMNS = ["h", "k", "l", "Amplitude", "weight", "Phase"] # Uses difference amplitudes
HKL_COLUMNS = ["h", "k", "l", "Amplitude", "Sigma"]
DARK_PHASE_COLUMNS = ["h", "k", "l", "Amplitude", "Filler", "Phase"]


def load_concentration_matrix(path):
    """
    Load the concentration/occupancy matrix.

    Returns
    -------
    labels : list[str]          timepoint labels, in file row order
    C : np.ndarray              shape (n_timepoints, n_species)
    species_names : list[str]   taken from the CSV header
    """
    df = pd.read_csv(path)
    label_col = df.columns[0]
    labels = df[label_col].astype(str).tolist()
    species_names = list(df.columns[1:])
    C = df[species_names].to_numpy(dtype=float)
    print(C)
    return labels, C, species_names


def load_phs(path):
    return pd.read_csv(path, sep=r"\s+", names=PHS_COLUMNS)


def load_hkl(path):
    return pd.read_csv(path, sep=r"\s+", names=HKL_COLUMNS)


def load_phase(path):
    return pd.read_csv(path, sep=r"\s+", names=DARK_PHASE_COLUMNS)


def adjust_ded_amplitudes_df(ded_df, dark_phase_df):
    """
    Adjust DED amplitudes if the phase differs from dark by 180 degrees.
    Uses vectorized pandas operations.
    """

    adjusted = ded_df.copy()

    # Keep only HKL and phase from dark
    dark_phase = dark_phase_df[["h", "k", "l", "Phase"]].copy()
    dark_phase = dark_phase.rename(columns={"Phase": "Phase_dark"})

    # Merge dark phases onto DED reflections
    adjusted = adjusted.merge(
        dark_phase,
        on=["h", "k", "l"],
        how="left"
    )

    # Calculate wrapped phase difference
    phase_ded = adjusted["Phase"] % 360
    phase_dark = adjusted["Phase_dark"] % 360

    delta = (phase_ded - phase_dark) % 360

    # Identify reflections requiring sign flip
    flip_mask = np.isclose(delta, 180, atol=1e-3)

    # Flip amplitudes
    adjusted.loc[flip_mask, "Amplitude"] *= -1

    # Rotate phases by 180° after flip
    adjusted.loc[flip_mask, "Phase"] = (
        phase_ded[flip_mask] + 180
    ) % 360

    # Normalize phases for all others
    adjusted.loc[~flip_mask, "Phase"] = phase_ded[~flip_mask]

    # Remove temporary column
    adjusted = adjusted.drop(columns=["Phase_dark"])

    print(f"Phase corrected reflections: {flip_mask.sum()}")

    return adjusted


def filter_common_reflections(data, data_sigma, min_occurrence):
    """
    Keep reflections present in at least `min_occurrence` timepoints.

    All output DataFrames are reindexed to the same HKL order.
    Missing reflections remain as NaN.
    """

    # Count HKL occurrences across timepoints
    hkl_counts = Counter()

    for df in data.values():
        hkl_counts.update(zip(df["h"], df["k"], df["l"]))

    valid_hkls = sorted(
        [hkl for hkl, count in hkl_counts.items()
         if count >= min_occurrence]
    )

    print(
        f"Keeping {len(valid_hkls)} reflections "
        f"present in at least {min_occurrence}/{len(data)} timepoints"
    )

    full_index = pd.MultiIndex.from_tuples(
        valid_hkls,
        names=["h", "k", "l"]
    )

    filtered_data = {}
    filtered_sigma = {}

    for label in data:
        # PHS
        df = data[label].set_index(["h", "k", "l"])
        filtered_data[label] = (
            df.reindex(full_index)
            .reset_index()
        )

        # HKL sigma
        df_sigma = data_sigma[label].set_index(["h", "k", "l"])
        filtered_sigma[label] = (
            df_sigma.reindex(full_index)
            .reset_index()
        )

    return filtered_data, filtered_sigma, full_index


def load_datasets(labels, phs_dir, phs_pattern, hkl_dir, hkl_pattern, dark_phase, min_occurrence):
    """
    Load per-timepoint difference amplitudes (.phs) and sigmas (.hkl).

    Only 'Amplitude' from the .phs files and 'Sigma' from the .hkl files are
    used numerically downstream. The .phs 'weight' column is recomputed
    later from propagated sigma, and the .phs h,k,l,Phase of the reference
    timepoint (see --reference-label) are reused for every output mode.
    """
    phs_dir = Path(phs_dir)
    hkl_dir = Path(hkl_dir)

    data, data_sigma = {}, {}

    for label in labels:
        data[label] = load_phs(phs_dir / phs_pattern.format(label=label))
        data[label] = adjust_ded_amplitudes_df(data[label], dark_phase)
        data_sigma[label] = load_hkl(hkl_dir / hkl_pattern.format(label=label))

    data, data_sigma, hkl_index = filter_common_reflections(
        data,
        data_sigma,
        min_occurrence=min_occurrence)

    return data, data_sigma, hkl_index


def build_matrices(labels, data, data_sigma, dark_sigma):
    """Stack per-timepoint amplitudes; propagate sigma by combining with dark sigma."""
    amplitudes_list = [data[label]["Amplitude"].to_numpy() for label in labels]
    data_matrix = np.vstack(amplitudes_list)

    sigma_list = [
        np.sqrt(data_sigma[label]["Sigma"].to_numpy() ** 2 + dark_sigma ** 2)
        for label in labels
    ]
    sigma_matrix = np.vstack(sigma_list)
    return data_matrix, sigma_matrix


def solve_kinetic_modes(C_matrix, data_matrix):
    """
    Solve C @ M ≈ D for M while handling NaNs in data_matrix.

    Parameters
    ----------
    C_matrix : np.ndarray
        Kinetic coefficient matrix (n_timepoints x n_species)
    data_matrix : np.ndarray
        Data matrix (n_timepoints x n_voxels)

    Returns
    -------
    M_basis : np.ndarray
        Kinetic modes (n_species x n_voxels)
    """

    n_timepoints, n_species = C_matrix.shape
    n_voxels = data_matrix.shape[1]

    M_basis = np.zeros((n_species, n_voxels))

    count = 0

    for hkl in range(n_voxels):

        y = data_matrix[:, hkl]

        # keep only timepoints with valid data
        idx = ~np.isnan(y)

        if np.sum(idx) < n_species:
            # not enough points to solve system
            count += 1
            continue

        C = C_matrix[idx, :]
        y = y[idx]

        try:
            M_basis[:, hkl] = np.linalg.solve(
                C.T @ C,
                C.T @ y
            )

        except Exception as e:
            print(
                "Valid timepoints:", np.sum(idx),
                "Rank:", np.linalg.matrix_rank(C),
                "Shape:", C.shape,
                "Error:", e
            )
            count += 1

    print(
        f"{count} reflections could not be decomposed because "
        f"the concentration matrix was rank-deficient."
    )

    return M_basis


def propagate_mode_sigma(species_weights, sigma_matrix):
    """
    Propagate reflection sigmas for one species, ignoring NaNs:
    sigma_j^2 = sum_t(C_tj^2 * sigma_t^2) / (sum_t C_tj)^2
    """

    valid = np.isfinite(sigma_matrix)

    w2 = (species_weights[:, None] ** 2) * valid
    weighted_sigma2 = np.where(valid, sigma_matrix**2, 0.0) * w2

    numerator = np.sum(weighted_sigma2, axis=0)
    denominator = np.sum(species_weights[:, None] * valid, axis=0) ** 2

    sigma_sq = np.full(sigma_matrix.shape[1], np.nan)
    mask = denominator > 0
    sigma_sq[mask] = numerator[mask] / denominator[mask]

    return np.sqrt(sigma_sq)


def voxel_weights_from_sigma(amplitude, sigma, eps=1e-12):
    """Down-weight reflections with large amplitude and/or large propagated sigma as proposed by Schmidt."""
    F2 = amplitude ** 2
    F2_mean = np.mean(F2)
    sigma2 = sigma ** 2
    sigma2_mean = np.mean(sigma2)
    w = 1.0 / (1.0 + F2 / (F2_mean + eps) + sigma2 / (sigma2_mean + eps))
    return w


def valid_voxel_mask(*arrays, cutoff):
    """True where every array is finite and at most `cutoff` in magnitude. Against numerical stability"""
    mask = np.ones_like(arrays[0], dtype=bool)
    for arr in arrays:
        mask &= ~np.isnan(arr) & (np.abs(arr) <= cutoff)
    return mask


def write_phs_file(df, output_path):
    with open(output_path, "w") as f:
        for _, row in df.iterrows():
            f.write(
                f"{int(row['h']):4d} {int(row['k']):4d} {int(row['l']):4d} "
                f"{row['Amplitude']:10.4f} {row['weight']:8.4f} {row['Phase']:8.3f}\n"
            )


def write_hkl_file(df, output_path):
    with open(output_path, "w") as f:
        for _, row in df.iterrows():
            f.write(
                f"{int(row['h']):4d} {int(row['k']):4d} {int(row['l']):4d} "
                f"{row['Amplitude']:10.4f} {row['Sigma']:8.4f}\n"
            )


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--concentrations", required=True,
                         help="CSV: timepoint label column + per-species occupancy columns")
    parser.add_argument("--phs-dir", required=True)
    parser.add_argument("--phs-pattern", default="{label}_diff.phs")
    parser.add_argument("--hkl-dir", required=True)
    parser.add_argument("--hkl-pattern", default="{label}_scaled.hkl")
    parser.add_argument("--dark-sigma", required=True, help=".hkl file with dark-state sigmas")
    parser.add_argument("--dark-phase", required=True, help=".hkl file with dark-state phases")
    parser.add_argument(
        "--amplitude-cutoff", type=float, default=1e6,
        help="Empirical cutoff: voxels with NaN or |amplitude|/sigma above this value "
             "are dropped from the output (guards against numerical instability in "
             "experimental data). Default: 1e6.",
    )
    parser.add_argument("--min-occurrence", type=int, required=True, help="Keep reflections present in at least this many timepoints")
    parser.add_argument("--out-prefix", default="state")
    args = parser.parse_args()

    labels, C, species_names = load_concentration_matrix(args.concentrations)
    if not 1 <= args.min_occurrence <= len(labels):
      parser.error(
          f"--min-occurrence must be between 1 and {len(labels)}"
      )
    print(f"Loaded concentration matrix: {C.shape[0]} timepoints x {C.shape[1]} species")
    print(f"Species: {species_names}")

    dark_phase = load_phase(args.dark_phase)

    data, data_sigma, hkl_index = load_datasets(
        labels, args.phs_dir, args.phs_pattern, args.hkl_dir, args.hkl_pattern, dark_phase, args.min_occurrence
    )

    dark_sigma_df = load_hkl(args.dark_sigma)

    dark_sigma = (dark_sigma_df.set_index(["h", "k", "l"]).reindex(hkl_index)["Sigma"].to_numpy())

    data_matrix, sigma_matrix = build_matrices(labels, data, data_sigma, dark_sigma)

    M_basis = solve_kinetic_modes(C, data_matrix)

    dark_phase_indexed = dark_phase.set_index(["h", "k", "l"])

    ref_df = pd.DataFrame(index=hkl_index).reset_index()

    ref_df["Phase"] = (dark_phase_indexed["Phase"].reindex(hkl_index).to_numpy())

    ref_df_sigma = data_sigma[labels[0]].copy()

    dark_ref_phase = ref_df["Phase"].to_numpy()

    for i, species in enumerate(species_names):
        mode_i = M_basis[i, :]

        # Sign correction: flip negative amplitudes and shift phase by 180 deg
        ref_phase = ref_df["Phase"].to_numpy()
        adjusted_amplitude = np.where(mode_i < 0, -mode_i, mode_i)
        adjusted_phase = np.where(mode_i < 0, (dark_ref_phase - 180), dark_ref_phase)

        sigma_propagated = propagate_mode_sigma(C[:, i], sigma_matrix)

        # --- .phs output: difference amplitude, phase, and a recomputed weight ---
        valid_phs = valid_voxel_mask(
            adjusted_amplitude, adjusted_phase, sigma_propagated, cutoff=args.amplitude_cutoff
        )
        weight = np.full_like(adjusted_amplitude, np.nan)
        weight[valid_phs] = voxel_weights_from_sigma(
            adjusted_amplitude[valid_phs], sigma_propagated[valid_phs]
        )

        df_phs = ref_df.copy()
        df_phs["Amplitude"] = adjusted_amplitude
        df_phs["Phase"] = adjusted_phase
        df_phs["weight"] = weight
        df_phs = df_phs.dropna(subset=["Amplitude", "Phase", "weight"])

        phs_path = f"{args.out_prefix}_{i + 1}_{species}.phs"
        write_phs_file(df_phs, phs_path)
        n_dropped_phs = len(adjusted_amplitude) - len(df_phs)
        print(
            f"Saved kinetic mode {i + 1} ({species}) to {phs_path} "
            f"({n_dropped_phs} voxels dropped by cutoff/NaN)"
        )

        # --- .hkl output: copied difference amplitude and propagated sigma ---
        raw_amplitude = np.abs(mode_i)
        valid_hkl = valid_voxel_mask(raw_amplitude, sigma_propagated, cutoff=args.amplitude_cutoff)

        df_hkl = ref_df_sigma.copy()
        df_hkl["Amplitude"] = raw_amplitude
        df_hkl["Sigma"] = sigma_propagated
        df_hkl = df_hkl[valid_hkl]

        hkl_path = f"{args.out_prefix}_{i + 1}_{species}_sigma.hkl"
        write_hkl_file(df_hkl, hkl_path)
        n_dropped_hkl = len(raw_amplitude) - len(df_hkl)
        print(
            f"Saved kinetic mode {i + 1} ({species}) sigma to {hkl_path} "
            f"({n_dropped_hkl} voxels dropped by cutoff/NaN)"
        )


if __name__ == "__main__":
    main()

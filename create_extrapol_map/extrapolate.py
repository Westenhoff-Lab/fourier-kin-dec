#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from collections import Counter

PHS_COLUMNS = ["h", "k", "l", "Amplitude_diff", "Weight", "Phase_diff"]
HKL_COLUMNS = ["h", "k", "l", "Amplitude", "Sigma"]
PHASE_COLUMNS = ["h", "k", "l", "Amplitude_ref", "Filler", "Phase_ref"]


def load_hkl(path):
    return pd.read_csv(path, sep=r"\s+", names=HKL_COLUMNS)


def load_phase_hkl(path):
    return pd.read_csv(path, sep=r"\s+", names=PHASE_COLUMNS)


def load_phs(path):
    return pd.read_csv(path, sep=r"\s+", names=PHS_COLUMNS)


def write_hkl(hkl, path):
    with open(path, "w") as f:
        for row in hkl.itertuples(index=False):
            f.write(
                f"{row.h:5d}{row.k:5d}{row.l:5d}"
                f"{row.Amplitude:14.4f}{row.Sigma:12.4f}{row.Phase:10.4f}\n"
            )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate weighted difference structure factors."
    )

    parser.add_argument(
        "--light",
        required=True,
        help="Light HKL file."
    )

    parser.add_argument(
        "--diff",
        required=True,
        help="Diff PHS file."
    )

    parser.add_argument(
        "--phase",
        required=True,
        help="Phase HKL file."
    )

    parser.add_argument(
        "--scale",
        type=float,
        required=True,
        help="Extrapolation factor S1."
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output difference structure factor file."
    )

    return parser.parse_args()

def main():

    args = parse_args()

    light_file = Path(args.light)
    phase_file = Path(args.phase)
    diff_file = Path(args.diff)
    output_file = Path(args.output)

    light = load_hkl(light_file)
    phase = load_phase_hkl(phase_file)
    diff = load_phs(diff_file)

    S1 = args.scale

    diff["Amplitude_diff"] = diff["Amplitude_diff"] * diff["Weight"] # We weight the difference amplitudes

    df = (
        phase
        .merge(
            diff[["h", "k", "l", "Amplitude_diff", "Phase_diff"]],
            on=["h", "k", "l"],
            how="left",
            sort=False
        )
        .merge(
            light[["h", "k", "l", "Sigma"]],
            on=["h", "k", "l"],
            how="left",
            sort=False
        )
    )
    df = df[
        (df["Amplitude_diff"].notna()) &
        (df["Sigma"].notna()) &
        (df["Amplitude_diff"] >= 0) &
        (df["Sigma"] >= 0)
    ]

    print(df)

    phase_difference = (
        (df["Phase_diff"] - df["Phase_ref"] + 180) % 360
    ) - 180

    same_direction = np.abs(phase_difference) < 20

    df["Fout"] = np.where(
        same_direction,
        df["Amplitude_ref"] + S1 * df["Amplitude_diff"],
        df["Amplitude_ref"] - S1 * df["Amplitude_diff"]
    )
    negative = df["Fout"] < 0

    df.loc[negative, "Fout"] *= -1

    df.loc[negative, "Phase_ref"] += 180

    df["Phase_ref"] = (
        (df["Phase_ref"] + 180) % 360
    ) - 180

    df["Sigma_out"] = np.sqrt(S1) * df["Sigma"]

    output = pd.DataFrame(
        {
            "h": df["h"],
            "k": df["k"],
            "l": df["l"],
            "Amplitude": df["Fout"],
            "Sigma": df["Sigma_out"],
            "Phase": df["Phase_ref"],
        }
    )

    write_hkl(output, output_file)

if __name__ == "__main__":
    main()

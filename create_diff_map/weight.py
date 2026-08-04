#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from collections import Counter

PHS_COLUMNS = ["h", "k", "l", "Amplitude", "Weight", "Phase"]
HKL_COLUMNS = ["h", "k", "l", "Amplitude", "Sigma"]
PHASE_COLUMNS = ["h", "k", "l", "Amplitude", "Sigma", "Phase"]


def load_hkl(path):
    return pd.read_csv(path, sep=r"\s+", names=HKL_COLUMNS)


def load_phase_hkl(path):
    return pd.read_csv(path, sep=r"\s+", names=PHASE_COLUMNS)


def write_phs(phs, path):
    with open(path, "w") as f:
        for row in phs.itertuples(index=False):
            f.write(
                f"{row.h:5d}{row.k:5d}{row.l:5d}"
                f"{row.Amplitude:10.4f}{row.Weight:10.4f}{row.Phase:10.4f}\n"
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
        "--dark",
        required=True,
        help="Dark HKL file."
    )

    parser.add_argument(
        "--phase",
        required=True,
        help="Phase HKL file."
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
    dark_file = Path(args.dark)
    phase_file = Path(args.phase)
    output_file = Path(args.output)

    light = load_hkl(light_file)
    dark = load_hkl(dark_file)
    phase = load_phase_hkl(phase_file)

    df = (
        light.rename(columns={"Amplitude": "F_light",
                              "Sigma": "SIGF_light"})
        .merge(
            dark.rename(columns={"Amplitude": "F_dark",
                                 "Sigma": "SIGF_dark"}),
            on=["h", "k", "l"],
            how="inner",
        )
        .merge(
            phase[["h", "k", "l", "Phase"]],
            on=["h", "k", "l"],
            how="inner",
        )
    )

    print(f"Matching reflections: {len(df)}")

    df["dF"] = df["F_light"] - df["F_dark"]
    df["sigma2"] = df["SIGF_light"]**2 + df["SIGF_dark"]**2

    DFSQM = np.mean(df["dF"]**2)
    S12SQM = np.mean(df["sigma2"])

    DFM = np.mean(df["dF"])
    ADFM = np.mean(np.abs(df["dF"]))
    S12M = np.mean(np.sqrt(df["sigma2"]))

    df["weight"] = 1.0 / (
        1.0
        + df["dF"]**2 / DFSQM
        + df["sigma2"] / S12SQM
    )

    WMEAN = df["weight"].mean()

    df["Phase_out"] = df["Phase"]
    neg = df["dF"] < 0
    df.loc[neg, "Phase_out"] = df.loc[neg, "Phase_out"] + 180.0
    df.loc[df["Phase_out"] > 180.0, "Phase_out"] -= 360.0

    df["dF_abs"] = np.abs(df["dF"])
    df["Amplitude"] = df["dF_abs"] / WMEAN

    print(df)

    phs = df[["h", "k", "l", "Amplitude", "weight", "Phase_out"]].copy()
    phs.columns = PHS_COLUMNS

    write_phs(phs, output_file)

if __name__ == "__main__":
    main()

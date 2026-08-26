"""Automate one KO run from input-deck generation through RMSE evaluation."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
from pathlib import Path

from postProcessKo import postProcessKo
from rmse import rmse
from writeKoDotIn import writeKoDotIn

WORKSHOP_DIR = Path(__file__).resolve().parent


def run_experiment(s1: float, yield_strength: float) -> float:
    """Generate KO input, run the solver, and return RMSE against the data."""
    compiler = shutil.which("gfortran")
    if compiler is None:
        message = "gfortran is required to run this workshop example."
        raise RuntimeError(message)

    required_files = ("KOv13.f", "TaData.txt")
    missing_files = [name for name in required_files if not (WORKSHOP_DIR / name).is_file()]
    if missing_files:
        message = f"Missing workshop input file(s): {', '.join(missing_files)}"
        raise FileNotFoundError(message)

    previous_directory = Path.cwd()
    try:
        os.chdir(WORKSHOP_DIR)
        input_path = writeKoDotIn(s1, yield_strength)
        print(f"Wrote {input_path.name} with s1={s1:.3f}, Y={yield_strength:.3e}")
        subprocess.run([compiler, "KOv13.f"], check=True)  # noqa: S603 - fixed local compiler/source
        subprocess.run(["./a.out"], check=True)
        postProcessKo()
        result = rmse("tv.txt", "TaData.txt")
    finally:
        os.chdir(previous_directory)

    return result


def parse_arguments() -> argparse.Namespace:
    """Return the two material parameters chosen for this run."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--s1", type=float, default=1.2, help="Mie-Gruneisen s1 parameter")
    parser.add_argument("--yield-strength", type=float, default=7.7e-3, help="Yield strength Y in Mbar")
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_arguments()
    rmse_value = run_experiment(arguments.s1, arguments.yield_strength)
    print(f"RMSE: {rmse_value:.6f} m/s")

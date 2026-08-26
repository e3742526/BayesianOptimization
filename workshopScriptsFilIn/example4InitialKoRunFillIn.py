"""Compile and run KO using the current workshop ``ko.in`` input deck."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

WORKSHOP_DIR = Path(__file__).resolve().parent


def run_ko() -> None:
    """Compile KO and run it from the directory containing its input deck."""
    compiler = shutil.which("gfortran")
    if compiler is None:
        message = "gfortran is required to run this workshop example."
        raise RuntimeError(message)

    required_files = ("KOv13.f", "ko.in")
    missing_files = [name for name in required_files if not (WORKSHOP_DIR / name).is_file()]
    if missing_files:
        message = f"Missing workshop input file(s): {', '.join(missing_files)}"
        raise FileNotFoundError(message)

    previous_directory = Path.cwd()
    try:
        os.chdir(WORKSHOP_DIR)
        subprocess.run([compiler, "KOv13.f"], check=True)  # noqa: S603 - fixed local compiler/source
        subprocess.run(["./a.out"], check=True)
    finally:
        os.chdir(previous_directory)


if __name__ == "__main__":
    run_ko()

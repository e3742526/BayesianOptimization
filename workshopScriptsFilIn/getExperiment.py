import subprocess
from pathlib import Path

from postProcessKo import postProcessKo
from rmse import rmse
from writeKoDotIn import writeKoDotIn

WORKSHOP_DIR = Path(__file__).resolve().parent


def getExperiment(s1,Y):

        writeKoDotIn(s1, Y)
        print(f"Generated ko.in with s1={s1} and Y={Y}")

        subprocess.run(
            ['gfortran', 'KOv13.f', '-o', 'a.out'],
            check=True,
            cwd=WORKSHOP_DIR,
        )
        print("Compilation successful!\n")

        subprocess.run(['./a.out'], check=True, cwd=WORKSHOP_DIR)

        postProcessKo(
            WORKSHOP_DIR / 'ko.dat',
            WORKSHOP_DIR / 'tv.txt',
        )
        print("Postprocessing successful!\n")

        simData = WORKSHOP_DIR / 'tv.txt'
        dataData = WORKSHOP_DIR / 'TaData.txt'

        rmseValue = rmse(simData, dataData)
        return rmseValue

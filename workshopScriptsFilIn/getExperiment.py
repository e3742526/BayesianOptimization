import subprocess
from postProcessKo import postProcessKo
from rmse import rmse
from writeKoDotIn import writeKoDotIn

def getExperiment(s1,Y):

        writeKoDotIn(s1, Y)
        print(f"Generated ko.in with s1={s1} and Y={Y}")

        subprocess.run(['gfortran', 'KOv13.f'], check=True)
        print("Compilation successful!\n")

        subprocess.run(['./a.out'], check=True)

        postProcessKo()
        print("Postprocessing successful!\n")

        simData = 'tv.txt'
        dataData = 'TaData.txt'

        rmseValue = rmse(simData, dataData)
        return rmseValue
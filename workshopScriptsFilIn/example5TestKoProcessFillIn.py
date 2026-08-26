from writeKoDotIn import writeKoDotIn
from postProcessKo import postProcessKo
from rmse import rmse
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import os

if __name__ == "__main__":
    import sys

    s1 = 1.2
    #s1Vec = np.linspace(0.5,3,10)

    plotValues = False
    
    Y = 7.7e-3
    if plotValues:
        YVec = np.linspace(0.005,0.015,11)
    else:
        YVec = np.array([Y])

    rmseOut = np.zeros(len(YVec)) 

    for i, Y in enumerate(YVec):

        ###################
        # Fill In Code Here
        ###################
        print(f"Generated ko.in with s1={s1} and Y={Y}")

        print("--- Running: gfortran KOv13.f ---")
        subprocess.run(['gfortran', 'KOv13.f'], check=True)
        print("Compilation successful!\n")

        # 4. Execute the compiled binary (./a.out)
        # Check if the file exists first to prevent execution errors
        subprocess.run(['./a.out'], check=True)

        ###################
        # Fill In Code Here
        ###################
        print("Postprocessing successful!\n")

        simData = 'tv.txt'
        dataData = 'TaData.txt'

        ###################
        # Fill In Code Here
        ###################
        rmseOut[i] = rmseValue
        print(f"RMSE between simulation and experimental data: {rmseValue:.6f}")


if plotValues:
    plt.rcParams.update({'font.size': 25})
    plt.plot(YVec, rmseOut, marker='o')
    plt.xlabel('Y')
    plt.ylabel('RMSE')
    plt.title('RMSE vs Y')
    plt.subplots_adjust(left=0.20, bottom=0.20)
    plt.show()
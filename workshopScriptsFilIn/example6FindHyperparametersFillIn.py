from writeKoDotIn import writeKoDotIn
from postProcessKo import postProcessKo
from rmse import rmse
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import os
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel
import matplotlib

if __name__ == "__main__":
    import sys

        ###################
        # Fill In Code Here
        ###################

    if runDoe:
        rmseOut = np.zeros(len(YVec)) 

    
        for i, Y in enumerate(YVec):

            writeKoDotIn(s1, Y)
            print(f"Generated ko.in with s1={s1} and Y={Y}")

            print("--- Running: gfortran KOv13.f ---")
            subprocess.run(['gfortran', 'KOv13.f'], check=True)
            print("Compilation successful!\n")

            # Execute the compiled binary (./a.out)
            # Check if the file exists first to prevent execution errors
            subprocess.run(['./a.out'], check=True)

            postProcessKo()
            print("Postprocessing successful!\n")

            simData = 'tv.txt'
            dataData = 'TaData.txt'

            rmseValue = rmse(simData, dataData)
            rmseOut[i] = rmseValue
            print(f"RMSE between simulation and experimental data: {rmseValue:.6f}")

        print(rmseOut)
    else:
         # I filled this in mannually
         rmse= np.array([36.65616066, 34.01236895, 49.47015776])

         X = YVec
         y = rmse
         matplotlib.rcParams.update({"font.size": 25})
         plt.figure(1)
         plt.plot(X, y ,'s', markersize=20)

         X = X.reshape(-1, 1)

         kernel =  ConstantKernel(10.0)*RBF(length_scale=0.1, length_scale_bounds='fixed')
         gp = GaussianProcessRegressor(
            kernel=kernel,
            alpha=1e-1,       # noise variance
            normalize_y=True)

         gp.fit(X, y)

         XPlot = np.linspace(min(X), max(X), 100).reshape(-1,1)
         yPred, yStd = gp.predict(XPlot, return_std=True)

         plt.fill_between(
            XPlot.flatten(),
            yPred.flatten()-1*yStd.flatten(),
            yPred.flatten()+1*yStd.flatten(),
            alpha=0.3
         )

        # GP mean
         plt.plot(
            XPlot,
            yPred,
            'b-',
            linewidth=2,
         )

         plt.xlabel('Yield Strength (Y), MBar')
         plt.ylabel('RMSE, m/s')
         plt.subplots_adjust(bottom=0.18, left=0.25)
         plt.show()


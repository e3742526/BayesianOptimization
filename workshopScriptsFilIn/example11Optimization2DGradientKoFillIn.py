from matplotlib.pylab import norm
import numpy as np
import matplotlib.pyplot as plt
import os
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel
import matplotlib
from scipy.stats import norm
import math
from scipy.stats import qmc
from scipy.optimize import minimize
import time
from getExperiment import getExperiment

start = time.perf_counter()

def runExperiment(X1,X2):
    s1 = X1
    Y = X2 
    if isinstance(Y, np.ndarray) and Y.ndim > 0: 
        output = np.zeros(len(Y))
        for i in range(len(Y)):
            output[i] = -1.0*getExperiment(s1[i], Y[i])
    else:
        output = -1.0*getExperiment(s1, Y)
    return output

def designExperiment(minVal, maxVal, numPoints):
    nearest = nearest = int(np.ceil(np.log2(numPoints)))
    sampler = qmc.Sobol(d=2, scramble=False)
    X = sampler.random_base2(m=nearest)   
    X = qmc.scale(X, [minVal[0], minVal[1]], [maxVal[0], maxVal[1]])
    Y = runExperiment(X[:,0],X[:,1])
    return X, Y

def acquisitionFunction(x, gp, acqFunType='UCB', yBest=None):
 
    # x should be shape (n_samples, n_features)
    yPred, yStd = gp.predict(x, return_std=True)
    if acqFunType == 'UCB':
        kappa = 0.1
        #kappa = 2.0
        acq = yPred + kappa * yStd
    elif acqFunType == 'PI':
        Z = (yPred - yBest) / yStd
        acq = norm.cdf(Z) 
    elif acqFunType == 'EI':  
        improvement = yPred - yBest
        Z = (yPred - yBest) / yStd
        acq = improvement * norm.cdf(Z) + yStd * norm.pdf(Z)    
    return acq

# ------------------------------------------------------------
# Acquisition function to maximize
# ------------------------------------------------------------
def objective(x):

    # acquisitionFunction expects an array of points
    X = np.asarray(x).reshape(1, -1)

    if acqFunType == 'PI' or acqFunType == 'EI':
        acq = acquisitionFunction(X, gp, acqFunType, yBest)
    else:
        acq = acquisitionFunction(X, gp, acqFunType)

    # L-BFGS-B minimizes, so return negative acquisition
    return -float(np.asarray(acq).ravel()[0])

def searchNextPoint(minVal, maxVal, gp, numPoints=16,
                    acqFunType='UCB', yBest=None):
    ###################
    # Fill In Code Here
    ###################

    return nextPoint

matplotlib.rcParams.update({"font.size": 25})

# s1 and 1
minVal = np.array([0.5,0.006])
maxVal = np.array([3,0.014])
numPoints = 2

#acqFunType = 'UCB'  # Upper Confidence Bound
#acqFunType = 'PI'   # Probability of Improvement
acqFunType = 'EI'   # Expected Improvement


XDOE, yDOE = designExperiment(minVal, maxVal, numPoints)

X = XDOE
# Ensure y is column vector (n_samples, 1) while X remains (n_samples, 2)
y = yDOE.reshape(-1, 1)

kernel =  ConstantKernel(1.0)*RBF(length_scale=[0.5, 0.001], length_scale_bounds='fixed')
gp = GaussianProcessRegressor(
    kernel=kernel,
    alpha=1e-6,       # noise variance
    normalize_y=True)

gp.fit(X, y)

nItr = 1

for i in range(nItr):
    if acqFunType == 'PI' or acqFunType == 'EI':
        yBest = np.max(y)
        XNew = searchNextPoint(minVal, maxVal, gp, 16, acqFunType,yBest)
    else:
        XNew = searchNextPoint(minVal, maxVal, gp, 16, acqFunType)
    # evaluate at the new point (XNew is a 2-element array)
    yNew = runExperiment(XNew[0], XNew[1])
    # append new sample (row) to X and new target to y
    X = np.vstack((X, XNew.reshape(1, -1)))
    y = np.vstack((y, np.array([[yNew]])))

    gp.fit(X, y)

    # create a 2D grid for plotting/prediction
    x1_plot = np.linspace(minVal[0], maxVal[0], 50)
    x2_plot = np.linspace(minVal[1], maxVal[1], 50)
    X1p, X2p = np.meshgrid(x1_plot, x2_plot)
    XPlot = np.column_stack((X1p.ravel(), X2p.ravel()))
    yPred, yStd = gp.predict(XPlot, return_std=True)
    yPred = np.asarray(yPred).reshape(X1p.shape)

bestPoint = X[np.argmax(y), :]

print('The optimum is ' + str(bestPoint))

# plot final
plt.figure(1)

x1p = np.linspace(minVal[0], maxVal[0], 50)
x2p = np.linspace(minVal[1] , maxVal[1], 50)
X1p, X2p = np.meshgrid(x1p, x2p)
Xp = np.column_stack((X1p.ravel(), X2p.ravel()))
yPred, yStd = gp.predict(Xp, return_std=True)
yPred = np.asarray(yPred).reshape(X1p.shape)
yStd = np.asarray(yStd).reshape(X1p.shape)
plt.contourf(X1p, X2p, yPred, cmap='jet', edgecolor='none')
plt.colorbar(label='y')
plt.plot(X[:,0],X[:,1],'ok')
plt.plot(XDOE[:,0],XDOE[:,1],'or')
plt.xlabel('$X_1$')
plt.ylabel('$X_2$')  
plt.subplots_adjust(bottom=0.25,  right=0.8, left=0.3)

elapsed = time.perf_counter() - start

print(f"Elapsed time: {elapsed:.3f} seconds")
print(f"Elapsed time: {elapsed/60:.2f} minutes")


plt.show()

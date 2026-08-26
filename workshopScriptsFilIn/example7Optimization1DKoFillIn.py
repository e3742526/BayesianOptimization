from matplotlib.pylab import norm
import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
from getExperiment import getExperiment
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel
import matplotlib
from scipy.stats import norm

WORKSHOP_DIR = Path(__file__).resolve().parent


def runExperiment(x):
    """Return negative RMSE for one or more yield-strength values."""
    s1 = 1.2
    y_values = np.atleast_1d(np.asarray(x, dtype=float))

    previous_directory = Path.cwd()
    try:
        os.chdir(WORKSHOP_DIR)
        output = np.array([-getExperiment(s1, y_value) for y_value in y_values])
    finally:
        os.chdir(previous_directory)

    return output.item() if np.ndim(x) == 0 else output

def designExperiment(minVal, maxVal, numPoints):
    rndX = np.linspace(minVal, maxVal, numPoints)
    yDOE = runExperiment(rndX)
    return rndX, yDOE

#Y is yield strength which is "x" in the slides

def acquisitionFunction(x, gp, acqFunType='UCB', yBest=None):
 
    yPred, yStd = gp.predict(x.reshape(-1,1), return_std=True)
    if acqFunType == 'UCB':
        kappa = 10.0
        acq = yPred + kappa * yStd
    elif acqFunType == 'PI':
        Z = (yPred - yBest) / yStd
        acq = norm.cdf(Z) 
    elif acqFunType == 'EI':  
        improvement = yPred - yBest
        Z = (yPred - yBest) / yStd
        acq = improvement * norm.cdf(Z) + yStd * norm.pdf(Z)    
    return acq

def searchNextPoint(minVal, maxVal, gp, numPoints=100, acqFunType='UCB', yBest=None):
    xSearch = np.linspace(minVal, maxVal, numPoints)
    if acqFunType == 'PI' or acqFunType == 'EI':
        acq = acquisitionFunction(xSearch, gp, acqFunType, yBest)
    else:
        acq = acquisitionFunction(xSearch, gp, acqFunType)
    nextPoint = xSearch[np.argmax(acq)]
    return nextPoint

minVal = 0.007
maxVal = 0.015
numPoints = 3

acqFunType = 'UCB'  # Upper Confidence Bound
#acqFunType = 'PI'   # Probability of Improvement
#acqFunType = 'EI'   # Expected Improvement

print('Designing initial experiment...')
XDOE, yDOE = designExperiment(minVal, maxVal, numPoints)
print('Experiment design complete.\n')

matplotlib.rcParams.update({"font.size": 25})
plt.figure(1)
plt.plot(XDOE, yDOE ,'s', markersize=20)

X = XDOE
y = yDOE

X = X.reshape(-1, 1)

kernel =  ConstantKernel(1.0)*RBF(length_scale=0.001, length_scale_bounds='fixed')
gp = GaussianProcessRegressor(
    kernel=kernel,
    alpha=1e-6,       # noise variance
    normalize_y=True)

print('Fitting Gaussian process...')
gp.fit(X, y)
print('Gaussian process fit complete.\n')

XPlot = np.linspace(minVal, maxVal, 100).reshape(-1,1)
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
plt.xlabel('X')
plt.ylabel('y')
plt.subplots_adjust(bottom=0.18, left=0.25)

nItr = 3

print(f"Running {nItr} iterations of Bayesian optimization...")
for i in range(nItr):
    if acqFunType == 'PI' or acqFunType == 'EI':
        yBest = np.max(y)
        XNew = searchNextPoint(minVal, maxVal, gp, 100, acqFunType,yBest)
    else:
        XNew = searchNextPoint(minVal, maxVal, gp, 100, acqFunType)
    yNew = runExperiment(XNew)
    X = np.append(X, XNew).reshape(-1,1)
    y = np.append(y, yNew)

    gp.fit(X, y)

    XPlot = np.linspace(minVal, maxVal, 100).reshape(-1,1)
    yPred, yStd = gp.predict(XPlot, return_std=True)

    print(f"Bayesian optimization {nItr} iterations completed.\n")

    XBest = X[np.argmax(y)]
    yBest = np.max(y)
    print(f"Best X: {XBest}, Best y: {yBest}\n")

plt.figure(2)
    
plt.plot(X,y, 'bs', markersize=10)
plt.xlim(0.0049,0.0151)
plt.ylim(-100,0.0)

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

plt.xlabel('X')
plt.ylabel('y')
plt.subplots_adjust(bottom=0.18, left=0.25)
    
print(f"Bayesian optimization completed.\n")

XBest = X[np.argmax(y)]
yBest = np.max(y)
print(f"Best X: {XBest}, Best y: {yBest}\n")

plt.show()

from matplotlib.pylab import norm
import numpy as np
import matplotlib.pyplot as plt
import os
try:
    import imageio
except Exception:
    imageio = None
try:
    from PIL import Image
    pil_available = True
except Exception:
    pil_available = False
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel
import matplotlib
from scipy.stats import norm

def runExperiment(x):
    sigma = 0.6
    frequency = 1.5
    output = (x) * np.exp(-0.5 * ((x-2.5) / sigma)**2) \
        * (1 + np.cos(2 * np.pi * frequency * (x-2.5))) / 2
    return output

def designExperiment(minVal, maxVal, numPoints):
    seed = 42;
    rng = np.random.default_rng(seed)
    #rndX = rng.uniform(minVal, maxVal, numPoints)
    rndX = np.linspace(minVal, maxVal, numPoints)
    yDOE = runExperiment(rndX)
    return rndX, yDOE

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

matplotlib.rcParams.update({"font.size": 25})
#bounds that it searches (minval to maxval))
minVal = 0.0
maxVal = 4.8
numPoints = 20# Design of experiments points

acqFunType = 'UCB'  # Upper Confidence Bound
#acqFunType = 'PI'   # Probability of Improvement
#acqFunType = 'EI'   # Expected Improvement

XDOE, yDOE = designExperiment(minVal, maxVal, numPoints)

# inital plot
plt.figure(1)
plt.plot(XDOE, yDOE ,'s', markersize=20)

X = XDOE
y = yDOE

X = X.reshape(-1, 1)

kernel =  ConstantKernel(10.0)*RBF(length_scale=0.25, length_scale_bounds='fixed')
gp = GaussianProcessRegressor(
    kernel=kernel,
    alpha=1e-4,       # noise variance
    normalize_y=True)

gp.fit(X, y)

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


nItr = 5
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


# final plot   
plt.figure(2)
    
plt.plot(X,y, 'bs', markersize=10)

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

# maximum point found from the experiment (not from the GP)
maxInd = np.argmax(y)
yMax = np.max(y)
plt.plot(X[maxInd],yMax,'rs', markersize=10)
print('The maximum is X = ' + str(X[maxInd]) + ' y = ' + str(yMax))

# plot the true solution (which you normally wouldn't know)
XTrue = np.linspace(minVal, maxVal, 100)
yTrue = runExperiment(XTrue)
plt.plot(XTrue, yTrue, 'g-', linewidth=2)

plt.xlabel('X')
plt.ylabel('y')
plt.subplots_adjust(bottom=0.18, left=0.25)
    
plt.show()

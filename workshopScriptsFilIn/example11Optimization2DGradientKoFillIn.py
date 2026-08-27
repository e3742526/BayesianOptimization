"""Workshop Example 11: fit two KO hydrocode parameters with Bayesian optimization.

Reference: AFRL Regional Net Workshop 2026, slides 109-110.  This is the
two-dimensional, expensive-experiment version of Example 10.  The input vector
is ``[s1, Y]``: ``s1`` ranges from 0.5 to 3 and ``Y`` from 0.006 to 0.014.
``getExperiment`` returns an RMSE-like mismatch to the Tantalum data, so this
file negates it and maximizes the resulting score.  Maximizing negative error
is exactly equivalent to minimizing the original error.

The marked ``searchNextPoint`` section is where the workshop combines a GP
acquisition function with bounded gradient optimization.  Every call that
reaches ``runExperiment`` may compile/run/post-process the KO model, so the
purpose of the surrogate is to decide carefully which single point to spend
that cost on next.
"""

# Kept from the workshop starter.  The later scipy.stats import provides the
# Normal CDF/PDF used by PI and EI, so it is the ``norm`` used in practice.
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

# Time the complete workflow, including initial KO evaluations and optimization.
start = time.perf_counter()

def runExperiment(X1,X2):
    # The optimizer represents a parameter set as [s1, Y].  Giving the
    # arguments descriptive local names makes the link to slide 109 explicit.
    s1 = X1
    Y = X2 
    # DOE evaluation passes vector columns, while the selected next point is a
    # scalar pair.  Supporting both forms lets one function serve both stages.
    if isinstance(Y, np.ndarray) and Y.ndim > 0: 
        output = np.zeros(len(Y))
        for i in range(len(Y)):
            # Negate RMSE because all acquisition functions and ``np.argmax``
            # in this file are written as maximizers.  A smaller true RMSE thus
            # becomes a larger, preferred objective value.
            output[i] = -1.0*getExperiment(s1[i], Y[i])
    else:
        output = -1.0*getExperiment(s1, Y)
    return output

def designExperiment(minVal, maxVal, numPoints):
    # ``random_base2`` needs a power-of-two sample count.  As in Example 10,
    # rounding up can generate more samples than the requested ``numPoints``.
    # The repeated assignment in the starter is harmless: both names receive
    # the same integer exponent used on the following line.
    nearest = nearest = int(np.ceil(np.log2(numPoints)))
    sampler = qmc.Sobol(d=2, scramble=False)
    X = sampler.random_base2(m=nearest)   
    # Scale unit-square Sobol samples to [0.5, 3] x [0.006, 0.014].  Preserve
    # the column order: X[:, 0] is s1 and X[:, 1] is Y throughout this script.
    X = qmc.scale(X, [minVal[0], minVal[1]], [maxVal[0], maxVal[1]])
    Y = runExperiment(X[:,0],X[:,1])
    return X, Y

def acquisitionFunction(x, gp, acqFunType='UCB', yBest=None):
    # x must be a table of candidate parameter pairs: (n_candidates, 2).
    # The GP returns a predicted negative-RMSE score and an uncertainty for each.
    yPred, yStd = gp.predict(x, return_std=True)
    if acqFunType == 'UCB':
        # UCB explicitly rewards uncertain regions.  Try the kappa values on
        # slide 110 to see the exploration/exploitation trade-off.
        kappa = 0.1
        #kappa = 2.0
        acq = yPred + kappa * yStd
    elif acqFunType == 'PI':
        # PI asks for the probability of beating the best observed score.
        Z = (yPred - yBest) / yStd
        acq = norm.cdf(Z) 
    elif acqFunType == 'EI':  
        # EI also accounts for the expected *size* of the improvement.
        improvement = yPred - yBest
        Z = (yPred - yBest) / yStd
        acq = improvement * norm.cdf(Z) + yStd * norm.pdf(Z)    
    return acq

def searchNextPoint(minVal, maxVal, gp, numPoints=16,
                    acqFunType='UCB', yBest=None):
    # Keep the optimizer callback local so it uses this call's GP, acquisition
    # type, and best observed value instead of unrelated module globals.
    def objective(x):
        # ``minimize`` supplies one vector; GP prediction expects a row batch.
        candidate = np.asarray(x).reshape(1, -1)
        acquisition = acquisitionFunction(
            candidate,
            gp,
            acqFunType,
            yBest,
        )

        # L-BFGS-B minimizes, so negate the acquisition we want to maximize.
        return -float(np.asarray(acquisition).ravel()[0])

    # Slide 110 calls for acquisition search plus gradient search.  Generate
    # bounded Sobol starts, run L-BFGS-B from each one, and retain the result
    # with the lowest objective (equivalently, the greatest acquisition).
    nearest = int(np.ceil(np.log2(numPoints)))
    sampler = qmc.Sobol(d=len(minVal), scramble=False)
    starts = sampler.random_base2(m=nearest)
    starts = qmc.scale(starts, minVal, maxVal)

    bounds = list(zip(minVal, maxVal))
    bestResult = None

    for startPoint in starts:
        result = minimize(
            objective,
            startPoint,
            method='L-BFGS-B',
            bounds=bounds,
        )

        if bestResult is None or result.fun < bestResult.fun:
            bestResult = result

    # ``numPoints`` always produces at least one Sobol start, so a best result
    # exists.  Return a plain two-element vector in the [s1, Y] order expected
    # by runExperiment and the append logic below.
    nextPoint = np.asarray(bestResult.x).reshape(-1)

    return nextPoint

matplotlib.rcParams.update({"font.size": 25})

# Bounds from slide 109.  The first coordinate is the KO s1 parameter and the
# second is yield strength Y (the comment above formerly said "1"; it means Y).
minVal = np.array([0.5,0.006])
maxVal = np.array([3,0.014])
numPoints = 2

#acqFunType = 'UCB'  # Upper Confidence Bound
#acqFunType = 'PI'   # Probability of Improvement
acqFunType = 'EI'   # Expected Improvement


# Initial Sobol samples provide the first expensive KO/RMSE observations.  The
# red points in the final plot identify this non-adaptive starting design.
XDOE, yDOE = designExperiment(minVal, maxVal, numPoints)

X = XDOE
# Ensure y is column vector (n_samples, 1) while X remains (n_samples, 2)
y = yDOE.reshape(-1, 1)

# Each input dimension has its own length scale because s1 and Y have very
# different numerical ranges.  A scalar length scale would imply the same
# distance scale in both units, which would be hard to interpret here.
kernel =  ConstantKernel(1.0)*RBF(length_scale=[0.5, 0.001], length_scale_bounds='fixed')
gp = GaussianProcessRegressor(
    kernel=kernel,
    alpha=1e-6,       # noise variance
    normalize_y=True)

# This GP maps parameter pairs to negative RMSE.  Higher predictions mean a
# better fit to the data because of the sign conversion in runExperiment.
gp.fit(X, y)

nItr = 1

for i in range(nItr):
    # PI/EI compare candidates to the best *measured* negative-RMSE score.
    # UCB only needs the GP mean and uncertainty.
    if acqFunType == 'PI' or acqFunType == 'EI':
        yBest = np.max(y)
        XNew = searchNextPoint(minVal, maxVal, gp, 16, acqFunType,yBest)
    else:
        XNew = searchNextPoint(minVal, maxVal, gp, 16, acqFunType)
    # The selected point is the one costly hydrocode evaluation in this loop.
    # All other acquisition/gradient operations use the inexpensive GP model.
    yNew = runExperiment(XNew[0], XNew[1])
    # append new sample (row) to X and new target to y
    X = np.vstack((X, XNew.reshape(1, -1)))
    y = np.vstack((y, np.array([[yNew]])))

    # Update the surrogate before the next search so it learns from the latest
    # KO result, even when the new result is worse than the current best.
    gp.fit(X, y)

    # Plot-only posterior grid; it does not add any new KO evaluations.
    x1_plot = np.linspace(minVal[0], maxVal[0], 50)
    x2_plot = np.linspace(minVal[1], maxVal[1], 50)
    X1p, X2p = np.meshgrid(x1_plot, x2_plot)
    XPlot = np.column_stack((X1p.ravel(), X2p.ravel()))
    yPred, yStd = gp.predict(XPlot, return_std=True)
    yPred = np.asarray(yPred).reshape(X1p.shape)

# Because y contains negative RMSE, the largest y corresponds to the *smallest*
# original RMSE.  This reports the best actually evaluated parameter pair.
bestPoint = X[np.argmax(y), :]

print('The optimum is ' + str(bestPoint))

# The contour is the final GP estimate of negative RMSE.  Black points are all
# KO evaluations; red points are just the initial DOE, as in Example 10.
fig, ax = plt.subplots(figsize=(8, 6), layout="constrained")

x1p = np.linspace(minVal[0], maxVal[0], 50)
x2p = np.linspace(minVal[1] , maxVal[1], 50)
X1p, X2p = np.meshgrid(x1p, x2p)
Xp = np.column_stack((X1p.ravel(), X2p.ravel()))
yPred, yStd = gp.predict(Xp, return_std=True)
yPred = np.asarray(yPred).reshape(X1p.shape)
yStd = np.asarray(yStd).reshape(X1p.shape)
posterior_contour = ax.contourf(X1p, X2p, yPred, cmap='jet', edgecolor='none')
fig.colorbar(posterior_contour, ax=ax, label='negative RMSE')
ax.plot(X[:,0],X[:,1],'ok', label="All KO evaluations")
ax.plot(XDOE[:,0],XDOE[:,1],'or', label="Initial DOE")
ax.set(title="Example 11: final GP estimate", xlabel='$s_1$', ylabel='$Y$')
ax.legend(loc="best")

elapsed = time.perf_counter() - start

print(f"Elapsed time: {elapsed:.3f} seconds")
print(f"Elapsed time: {elapsed/60:.2f} minutes")


plt.show()

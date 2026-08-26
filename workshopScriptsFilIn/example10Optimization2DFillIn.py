"""Workshop Example 10: Bayesian optimization of a 2-D teaching function.

Reference: AFRL Regional Net Workshop 2026, slides 93-94.  A Sobol DOE first
creates an initial GP model.  Bayesian optimization then repeatedly scores
candidate inputs with an acquisition function, evaluates the most promising
candidate, and refits the GP.  The exercise is to implement the candidate
search in ``searchNextPoint`` and compare UCB, PI, and EI.

This script *maximizes* ``runExperiment``.  That convention matters for all
three acquisition functions and for the final ``np.argmax``.  Example 11 will
reuse the same loop but negate an RMSE so that minimizing error becomes a
maximization problem.
"""

# ``scipy.stats.norm`` below is the Normal distribution used for PI/EI.  This
# earlier workshop import is retained unchanged; the later SciPy import shadows
# its name before any acquisition calculation is made.
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

import time

# Start timing before the DOE so the reported time covers the whole optimization
# workflow.  For this analytic example it should be short; the same structure
# becomes significant when each experiment invokes the KO hydrocode.
start = time.perf_counter()

def runExperiment(X1,X2):
    # A known smooth objective lets us test optimization logic cheaply.  NumPy
    # accepts scalar inputs (one new point) and arrays (the whole initial DOE).
    amp = 0.1
    Y = amp*np.sin(X1)  + amp*np.cos(X2)  
    return Y

def designExperiment(minVal, maxVal, numPoints):
    # Sobol's base-2 API needs 2**m samples.  Rounding *up* means the requested
    # 12 points on slide 93 becomes 16 actual DOE points (ceil(log2(12)) = 4).
    # This is a useful detail to verify whenever you change ``numPoints``.
    nearest =  int(np.ceil(np.log2(numPoints)))
    sampler = qmc.Sobol(d=2, scramble=False)
    X = sampler.random_base2(m=nearest)   
    # Map Sobol's unit-square coordinates to the physical search rectangle.
    # Each row is [X_1, X_2], matching the feature-column convention of Ex. 8-9.
    X = qmc.scale(X, [minVal[0], minVal[1]], [maxVal[0], maxVal[1]])
    Y = runExperiment(X[:,0],X[:,1])
    return X, Y

def acquisitionFunction(x, gp, acqFunType='UCB', yBest=None):
    # The GP gives a posterior mean (exploitation signal) and standard deviation
    # (exploration signal) for every candidate row in x.
    yPred, yStd = gp.predict(x, return_std=True)
    if acqFunType == 'UCB':
        # Upper Confidence Bound maximizes mean + kappa*uncertainty.  Larger
        # kappa values explore more; slide 94 suggests trying 0.1, 1, and 10.
        kappa = 1.0
        acq = yPred + kappa * yStd
    elif acqFunType == 'PI':
        # Probability of Improvement asks only "how likely is a result better
        # than yBest?"  It can be greedy because it ignores *how much* better.
        Z = (yPred - yBest) / yStd
        acq = norm.cdf(Z) 
    elif acqFunType == 'EI':  
        # Expected Improvement weights both probability and magnitude of a
        # potential improvement.  ``norm.cdf``/``norm.pdf`` come from SciPy.
        improvement = yPred - yBest
        Z = (yPred - yBest) / yStd
        acq = improvement * norm.cdf(Z) + yStd * norm.pdf(Z)    
    return acq

def searchNextPoint(minVal, maxVal, gp, numPoints=100, acqFunType='UCB', yBest=None):
    ###################
    # Fill In Code Here
    # Slide 93 specifies 100 Sobol candidates across [-pi, pi]^2.  A direct
    # implementation should: (1) generate candidate rows in that rectangle,
    # (2) evaluate ``acquisitionFunction`` for every row, and (3) return the
    # *single* row with the largest score.  The acquisition is cheap because it
    # queries the GP; do not call ``runExperiment`` for every candidate.
    #
    # ``nextPoint`` must be a two-element array, not shape (1, 2), because the
    # loop below indexes it as XNew[0] and XNew[1] before reshaping it for vstack.
    ###################
    return nextPoint

matplotlib.rcParams.update({"font.size": 25})

# A two-element bound vector makes every later operation visibly two-dimensional:
# minVal[0]/maxVal[0] belong to X_1 and minVal[1]/maxVal[1] belong to X_2.
minVal = np.array([-np.pi,-np.pi])
maxVal = np.array([np.pi,np.pi])
numPoints = 12

acqFunType = 'UCB'  # Upper Confidence Bound
#acqFunType = 'PI'   # Probability of Improvement
#acqFunType = 'EI'   # Expected Improvement


# The DOE establishes an initial surrogate before acquisition-driven decisions
# begin.  Red points in the final figure are these initial samples.
XDOE, yDOE = designExperiment(minVal, maxVal, numPoints)

X = XDOE
# Scikit-learn accepts a 1-D target too, but the column form makes appending
# each scalar yNew unambiguous and matches the earlier workshop examples.
y = yDOE.reshape(-1, 1)

kernel =  ConstantKernel(1.0)*RBF(length_scale=1.0, length_scale_bounds='fixed')
gp = GaussianProcessRegressor(
    kernel=kernel,
    alpha=1e-6,       # noise variance
    normalize_y=True)

# Fit the initial GP.  Each later iteration refits it after adding exactly one
# true experiment evaluation, which is the update loop in Bayesian optimization.
gp.fit(X, y)

nItr = 10

for i in range(nItr):
    # PI and EI need the best observed objective so far; UCB does not.  This
    # uses actual evaluations in y, never the GP's predicted maximum.
    if acqFunType == 'PI' or acqFunType == 'EI':
        yBest = np.max(y)
        XNew = searchNextPoint(minVal, maxVal, gp, 100, acqFunType,yBest)
    else:
        XNew = searchNextPoint(minVal, maxVal, gp, 100, acqFunType)
    # This is the expensive step in a real application.  The candidate search
    # above may inspect 100 GP estimates, but only this one location is sent to
    # the actual experiment in each iteration.
    yNew = runExperiment(XNew[0], XNew[1])
    # append new sample (row) to X and new target to y
    X = np.vstack((X, XNew.reshape(1, -1)))
    y = np.vstack((y, np.array([[yNew]])))

    # Refit after augmenting the training set so the next acquisition search
    # incorporates what the latest experiment taught us.
    gp.fit(X, y)

    # These predictions are useful for debugging/visualization.  They do not
    # select the next point; selection is controlled by searchNextPoint above.
    x1_plot = np.linspace(minVal[0], maxVal[0], 50)
    x2_plot = np.linspace(minVal[1], maxVal[1], 50)
    X1p, X2p = np.meshgrid(x1_plot, x2_plot)
    XPlot = np.column_stack((X1p.ravel(), X2p.ravel()))
    yPred, yStd = gp.predict(XPlot, return_std=True)
    yPred = np.asarray(yPred).reshape(X1p.shape)

# ``np.argmax`` finds the largest *observed* y, not necessarily the maximum of
# the smooth GP contour.  That distinction is important: this is the best
# validated experiment the algorithm has actually run.
bestPoint = X[np.argmax(y), :]

print('The optimum is ' + str(bestPoint))

# In the final figure, black points are all evaluated locations and red points
# identify the initial DOE subset.  The colored surface is the final GP mean.
fig, ax = plt.subplots(figsize=(8, 6), layout="constrained")
x1p = np.linspace(-np.pi, np.pi, 50)
x2p = np.linspace(-np.pi, np.pi, 50)
X1p, X2p = np.meshgrid(x1p, x2p)
Xp = np.column_stack((X1p.ravel(), X2p.ravel()))
yPred, yStd = gp.predict(Xp, return_std=True)
yPred = np.asarray(yPred).reshape(X1p.shape)
yStd = np.asarray(yStd).reshape(X1p.shape)
posterior_contour = ax.contourf(X1p, X2p, yPred, cmap='jet', edgecolor='none')
fig.colorbar(posterior_contour, ax=ax, label='y')
ax.plot(X[:,0],X[:,1],'ok', label="All evaluations")
ax.plot(XDOE[:,0],XDOE[:,1],'or', label="Initial DOE")
ax.set(title="Example 10: final GP posterior mean", xlabel='$X_1$', ylabel='$X_2$')
ax.legend(loc="best")
   
elapsed = time.perf_counter() - start

print(f"Elapsed time: {elapsed:.3f} seconds")
print(f"Elapsed time: {elapsed/60:.2f} minutes")
plt.show()

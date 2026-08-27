"""Workshop Example 9: compare two-dimensional designs of experiments.

Reference: AFRL Regional Net Workshop 2026, slides 90-91.  Example 8 gave the
GP every point on a dense grid.  This example deliberately keeps only a small
set of experiment locations, then asks how well those locations span the
two-dimensional domain and support a surrogate prediction.

The marked section is the experiment-design exercise.  It should generate
``X`` with shape ``(2**nDoeExp, 2)`` in the physical domain [-pi, pi]^2.  The
remaining code makes the same ``X``/``Y`` shape conversions introduced in
Example 8 and visualizes the result.
"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel
from scipy.stats import qmc
from scipy.interpolate import griddata

def runExperiment(X1,X2):
    # This analytic function stands in for a costly experiment.  NumPy applies
    # sin/cos elementwise, so X1 and X2 can be either scalars or vectors of DOE
    # locations.  That vectorization evaluates every DOE point in one call.
    amp = 0.1
    Y = amp*np.sin(X1)  + amp*np.cos(X2)  
    return Y

# Sobol's ``random_base2`` generator is designed for powers of two.  The slide
# asks you to increase this exponent and judge when the DOE becomes useful:
# nDoeExp=2 gives 4 points, nDoeExp=3 gives 8, and so on.
nDoeExp = 9 #(nPoints = 2^nDoeExp)

### different types of experiments below
doeType = 'sobol'
#doeType = 'LHC' # lattice hypercube

#doeType = 'random'

###################
# Fill In Code Here
if doeType == 'sobol':
    sampler = qmc.Sobol(d=2, scramble=False) 
    # Set to False, will give the same points every time, but scramble=True is more typical.
    X = sampler.random_base2(m=nDoeExp)
    X = qmc.scale(X, [-np.pi, -np.pi], [np.pi, np.pi])
elif doeType == 'LHC':
    sampler = qmc.LatinHypercube(d=2)
    X = sampler.random(n=2**nDoeExp)
    X = qmc.scale(X, [-np.pi, -np.pi], [np.pi, np.pi])
elif doeType == 'random':
    X = np.random.uniform(low=-np.pi, high=np.pi, size=(2**nDoeExp,2))
else:
    print('doeType not recognized') 

# sobol, LHC, random
###################
# Implement one of the DOE choices from slide 91 here:
#
# * Sobol samples are low-discrepancy points.  ``random_base2(m=nDoeExp)``
#   returns exactly 2**nDoeExp points in the unit square [0, 1)^2; use
#   ``qmc.scale`` to map both columns to [-pi, pi].
# * Latin hypercube sampling (LHC) stratifies each one-dimensional projection.
#   It does not require a power-of-two count, but use the same count here for a
#   fair visual comparison.
# * Random sampling is the baseline: it may leave clusters and holes simply by
#   chance.  The slide's question is which approach covers the space best.
#
# Set the final result to ``X`` with two columns.  Column 0 is X_1 and column 1
# is X_2; all downstream code relies on that convention.

# Split the two input-feature columns only because the analytic function and
# plotting APIs use separate coordinate arrays.  Keep X itself for GP fitting.
X1 = X[:,0]
X2 = X[:,1]
Y = runExperiment(X1,X2)

# This dense grid is for visualization, not for running new experiments.
# Interpolation below estimates a smooth-looking surface from the sparse DOE
# values, whereas the GP later supplies a probabilistic surrogate.
x1_grid = np.linspace(-np.pi, np.pi, 100)
x2_grid = np.linspace(-np.pi, np.pi, 100)

X1_grid, X2_grid = np.meshgrid(x1_grid, x2_grid)

# ``griddata`` is a deterministic visualization technique: it fills a regular
# grid from scattered data, but it does not quantify uncertainty.  Comparing
# this plot with the GP plot is a useful way to see why DOE coverage matters.
Y_grid = griddata(
    (X1, X2),
    Y,
    (X1_grid, X2_grid),
    method='cubic'
)

# Black dots identify the only locations actually evaluated.  Empty/NaN areas
# near a sparse convex hull are a normal interpolation limitation, not GP data.
matplotlib.rcParams.update({"font.size": 18})
# The interpolation and GP prediction answer complementary questions, so place
# them side by side in one window.  Constrained layout keeps both colorbars and
# labels visible without manually adjusting margins for separate figures.
fig, axes = plt.subplots(
    1, 2, figsize=(13, 5.5), layout="constrained", sharex=True, sharey=True
)
fig.suptitle(f"Example 9: {doeType} design of experiments")

interpolation_contour = axes[0].contourf(X1_grid, X2_grid, Y_grid, levels=30, cmap='jet')
fig.colorbar(interpolation_contour, ax=axes[0], label='y')
axes[0].scatter(X1, X2, color='black', label="DOE samples")
axes[0].set(title="Scattered-data interpolation", xlabel='$X_1$', ylabel='$X_2$')

# Build the sample-by-feature table expected by scikit-learn.  Ravel is harmless
# here because X1/X2 are already one-dimensional, and makes the convention
# explicit for readers carrying it forward from the gridded Example 8.
X = np.column_stack((X1.ravel(), X2.ravel()))
Y = Y.reshape(-1, 1)
print('X is:', X.shape)
print('Y is:', Y.shape)
# A length scale of one is a workshop choice, not a universal optimum.  It
# means changes are expected on an O(1) distance in input space.  Slide 91
# encourages changing the DOE while holding this modeling choice fixed, so
# differences in the plots can be attributed mainly to experiment design.
kernel =  ConstantKernel(1.0)*RBF(length_scale=1, length_scale_bounds='fixed')

gp = GaussianProcessRegressor(
    kernel=kernel,
    alpha=1e-6,       # noise variance
    normalize_y=True)

# Fit only to the DOE points; unlike Example 8, the GP never sees the dense
# 100 x 100 visualization grid or the 50 x 50 prediction grid.
gp.fit(X, Y)

# Make a separate 50 x 50 regular grid on which to inspect the posterior mean.
# ``Xp`` has 2,500 two-feature rows; reshape returns the estimates to plot form.
x1p = np.linspace(-np.pi, np.pi, 50)
x2p = np.linspace(-np.pi, np.pi, 50)
X1p, X2p = np.meshgrid(x1p, x2p)
Xp = np.column_stack((X1p.ravel(), X2p.ravel()))

# The mean is plotted below.  ``yStd`` is retained so you can easily add the
# uncertainty plot from Example 8 and relate uncertainty to DOE coverage.
yPred, yStd = gp.predict(Xp, return_std=True)

# Reshape predictions to match the grid shape (X1, X2)
yPred = np.asarray(yPred).reshape(X1p.shape)
yStd = np.asarray(yStd).reshape(X1p.shape)

prediction_contour = axes[1].contourf(X1p, X2p, yPred, cmap='jet', edgecolor='none')
fig.colorbar(prediction_contour, ax=axes[1], label='y')
axes[1].scatter(X1, X2, color='black', label="DOE samples")
axes[1].set(title="GP posterior mean", xlabel='$X_1$', ylabel='$X_2$')

plt.show()

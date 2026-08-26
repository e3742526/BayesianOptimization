"""Workshop Example 8: represent a two-dimensional function for a GP.

Reference: AFRL Regional Net Workshop 2026, slides 87-88 ("Example 8: 2D
Example" and "Let's Code/Try It").  The central lesson is that a plotting
grid and a machine-learning data table describe the same locations in two
different shapes:

* ``X1`` and ``X2`` are two 2-D arrays used to draw a contour plot.
* ``X`` is a 2-D table of samples with shape ``(n_samples, 2)``; this is the
  shape ``GaussianProcessRegressor`` expects.
* ``Y`` is the corresponding target table with shape ``(n_samples, 1)``.
"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel

minVal = -1.0*np.pi
maxVal = 1.0*np.pi

# The workshop uses the square domain -pi < X_1, X_2 < pi.  Both coordinate
# directions use the same bounds, but they remain distinct input features once
# the grid is converted to ``X``.
#
# Slide 87 asks for a 50 x 50 grid.  ``linspace`` creates each axis, while
# ``meshgrid`` makes arrays of all coordinate pairs.  At index [row, column],
# ``(X1[row, column], X2[row, column])`` is one location in the 2-D domain.
x1 = np.linspace(minVal, maxVal, 50)
x2 = np.linspace(minVal, maxVal, 50)
X1, X2 = np.meshgrid(x1, x2) #combos of all the x1 and x2 values
# The GP does not receive a mesh.  It receives one row per experiment, with
# columns [X_1, X_2].  Raveling both grids in the same order preserves the
# coordinate pairs, producing the slide's expected shape: (2500, 2).
X = np.column_stack((X1.ravel(), X2.ravel()))

# This inexpensive analytic function is the known "truth" for the exercise.
# Later KO examples replace it conceptually with a costly simulation.  The
# amplitude controls the output scale, not the location of the maximum.
amp = 0.1
Y = amp*np.sin(X1)  + amp*np.cos(X2)  

# This is a filled contour plot (a top-down view), not a 3-D perspective plot.
# It visualizes the known function before comparing it with the GP prediction.
matplotlib.rcParams.update({"font.size": 18})
# Keep the three related views in one window.  ``layout="constrained"`` lets
# Matplotlib reserve room for labels and each colorbar without manual margin
# tuning or three independently stacking figure windows.
fig, axes = plt.subplots(1, 3, figsize=(18, 5.5), layout="constrained")
fig.suptitle("Example 8: true function, GP mean, and GP uncertainty")

truth_contour = axes[0].contourf(X1, X2, Y,cmap='jet', edgecolor='none')
fig.colorbar(truth_contour, ax=axes[0], label='y')
axes[0].set(title="Known function", xlabel='$X_1$', ylabel='$X_2$')

###################
# Fill In Code Here
# x1, x2, X1,X2
###################
# There is no second grid construction needed after the first section is
# complete.  This retained workshop marker is a good checkpoint: verify that
# X has one row per Y value and two input-feature columns before fitting.
Y = Y.reshape(-1, 1)
# ``-1`` asks NumPy to infer the row count.  Here it turns 2,500 values into a
# one-column target table.  Each Y row must correspond to the X row with the
# same index, which the shared ravel order above guarantees.
print('X is:', X.shape)
print('Y is:', Y.shape)

# The kernel encodes the belief that nearby inputs have correlated outputs.  The
# RBF length scale says what "nearby" means in X_1/X_2 units.  Slide 88 asks
# you to try length scales and inspect uncertainty; it is fixed here so the
# exercise changes a visible modeling choice rather than learning it silently.
kernel =  ConstantKernel(1.0)*RBF(length_scale=0.01, length_scale_bounds='fixed')
gp = GaussianProcessRegressor(
    kernel=kernel,
    alpha=1e-6,       # noise variance
    normalize_y=True)

# Training stores the relationship between all 2,500 locations and their known
# values.  This dense reference case is deliberately easier than Example 9,
# where the GP receives only a sparse design of experiments.
gp.fit(X, Y)

# Slide 87 then asks for predictions on a new, coarser 10 x 10 grid.  Maintain
# both representations: X1p/X2p for plotting and Xp as 100 GP input rows.
x1p = np.linspace(minVal, maxVal, 10)
x2p = np.linspace(minVal, maxVal, 10)
X1p, X2p = np.meshgrid(x1p, x2p)

Xp = np.column_stack((X1p.ravel(), X2p.ravel()))

# ``yPred`` is the posterior mean (the GP's best estimate).  ``yStd`` is the
# posterior standard deviation: uncertainty in the model, not noise measured
# from repeated experiments.  High uncertainty identifies where a new sample
# could be informative.
yPred, yStd = gp.predict(Xp, return_std=True)

# GP output is returned in table order.  Reshape it back to the 10 x 10 grid so
# contourf can place each prediction at its original coordinate.
yPred = np.asarray(yPred).reshape(X1p.shape)
yStd = np.asarray(yStd).reshape(X1p.shape)

prediction_contour = axes[1].contourf(X1p, X2p, yPred, cmap='jet', edgecolor='none')
fig.colorbar(prediction_contour, ax=axes[1], label='y')
axes[1].scatter(X1p, X2p, color='black', s=12, label="Prediction grid")
axes[1].set(title="GP posterior mean", xlabel='$X_1$', ylabel='$X_2$')

# This addresses the slide-88 uncertainty question.  With the very short 0.01
# length scale, observations influence only tiny neighborhoods, so large areas
# of the domain can remain uncertain despite the regular training grid.
uncertainty_contour = axes[2].contourf(X1p, X2p, yStd, cmap='jet', edgecolor='none')
fig.colorbar(uncertainty_contour, ax=axes[2], label=r'$\sigma$')
axes[2].set(title="GP uncertainty", xlabel='$X_1$', ylabel='$X_2$')


plt.show()

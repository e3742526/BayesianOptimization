import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from sklearn.neural_network import MLPRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel

# Load the values printed on slide 50.
data_path = Path(__file__).with_name('example2Slide50Data.txt')
data = np.loadtxt(data_path, delimiter=',')
X = data[:,0]
y = data[:,1]

# make column data
X = X.reshape(-1, 1)


# Kernel and GP setup
kernel = ConstantKernel(1) * RBF(1)
gp = GaussianProcessRegressor(
    kernel=kernel,
    alpha=0.2,  # noise variance
    normalize_y=True,
)

# Train GP
gp.fit(X, y)

# Extract fitted RBF length scale
fitted_kernel = gp.kernel_
if hasattr(fitted_kernel, 'length_scale'):
    fitted_length_scale = fitted_kernel.length_scale
elif hasattr(fitted_kernel, 'k2') and hasattr(fitted_kernel.k2, 'length_scale'):
    fitted_length_scale = fitted_kernel.k2.length_scale
else:
    fitted_length_scale = None
print(f"Fitted RBF length scale: {fitted_length_scale}")

# Points for interpolation
Xstar = np.linspace(X.min(), X.max(), 200).reshape(-1, 1)

# Mean and standard deviation
mean, std = gp.predict( Xstar, return_std=True)

# Format figure
matplotlib.rcParams.update({"font.size": 25})
plt.figure(1)

# Plot data
plt.plot(X, y, 'rs', markersize=16)

# GP Standard Deviation Bounds
plt.fill_between(
    Xstar.flatten(),
    mean.flatten()-std.flatten(),
    mean.flatten()+std.flatten(),
    alpha=0.3
)

# GP mean
plt.plot(
    Xstar,
    mean,
    'b-',
    linewidth=2,
)

# label figure
plt.xlabel('Parameter')
plt.ylabel('Blast Resistance')
plt.subplots_adjust(bottom=0.18, left=0.25)

plt.show()

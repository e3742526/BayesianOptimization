import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel
from scipy.stats import qmc
from scipy.interpolate import griddata

def runExperiment(X1,X2):
    amp = 0.1
    Y = amp*np.sin(X1)  + amp*np.cos(X2)  
    return Y

nDoeExp = 2 #(nPoints = 2^nDoeExp)

doeType = 'sobol'
#doeType = 'LHC' # lattice hypercube
#doeType = 'random'

###################
# Fill In Code Here
# sobol, LHC, random
###################

X1 = X[:,0]
X2 = X[:,1]
Y = runExperiment(X1,X2)

x1_grid = np.linspace(-np.pi, np.pi, 100)
x2_grid = np.linspace(-np.pi, np.pi, 100)

X1_grid, X2_grid = np.meshgrid(x1_grid, x2_grid)

# Interpolate Y onto the grid
Y_grid = griddata(
    (X1, X2),
    Y,
    (X1_grid, X2_grid),
    method='cubic'
)

# Plot interpolated contour
matplotlib.rcParams.update({"font.size": 18})
plt.figure(1)
plt.contourf(X1_grid, X2_grid, Y_grid, levels=30, cmap='jet')
plt.colorbar(label='y')
plt.scatter(X1, X2, color='black')
plt.xlabel('$X_1$')
plt.ylabel('$X_2$')  
plt.subplots_adjust(bottom=0.18,  right=0.9) 

X = np.column_stack((X1.ravel(), X2.ravel()))
Y = Y.reshape(-1, 1)
print('X is:', X.shape)
print('Y is:', Y.shape)
kernel =  ConstantKernel(1.0)*RBF(length_scale=1, length_scale_bounds='fixed')

gp = GaussianProcessRegressor(
    kernel=kernel,
    alpha=1e-6,       # noise variance
    normalize_y=True)

gp.fit(X, Y)

# inputs for prediction
x1p = np.linspace(-np.pi, np.pi, 50)
x2p = np.linspace(-np.pi, np.pi, 50)
X1p, X2p = np.meshgrid(x1p, x2p)
Xp = np.column_stack((X1p.ravel(), X2p.ravel()))

# prediction
yPred, yStd = gp.predict(Xp, return_std=True)

# Reshape predictions to match the grid shape (X1, X2)
yPred = np.asarray(yPred).reshape(X1p.shape)
yStd = np.asarray(yStd).reshape(X1p.shape)

plt.figure(2)
plt.contourf(X1p, X2p, yPred, cmap='jet', edgecolor='none')
plt.colorbar(label='y')
plt.scatter(X1, X2, color='black')
plt.xlabel('$X_1$')
plt.ylabel('$X_2$')  
plt.subplots_adjust(bottom=0.18,  right=0.9)

plt.show()
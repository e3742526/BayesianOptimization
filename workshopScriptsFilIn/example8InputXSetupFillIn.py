import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel

minVal = -1.0*np.pi
maxVal = 1.0*np.pi

# Create grid coordinates (similar to MATLAB's default range)
###################
# Fill In Code Here
# x1, x2, X1,X2
###################

amp = 0.1
Y = amp*np.sin(X1)  + amp*np.cos(X2)  

# Plot the 3D surface
matplotlib.rcParams.update({"font.size": 18})
plt.figure(1)
plt.contourf(X1, X2, Y,cmap='jet', edgecolor='none')
plt.colorbar(label='y')
plt.xlabel('$X_1$')
plt.ylabel('$X_2$')  
plt.subplots_adjust(bottom=0.18,  right=0.8,left=0.25) 

###################
# Fill In Code Here
# x1, x2, X1,X2
###################
Y = Y.reshape(-1, 1)
print('X is:', X.shape)
print('Y is:', Y.shape)

kernel =  ConstantKernel(1.0)*RBF(length_scale=0.01, length_scale_bounds='fixed')
gp = GaussianProcessRegressor(
    kernel=kernel,
    alpha=1e-6,       # noise variance
    normalize_y=True)

gp.fit(X, Y)

# Predict the result at 10 x 10 new points
x1p = np.linspace(minVal, maxVal, 10)
x2p = np.linspace(minVal, maxVal, 10)
X1p, X2p = np.meshgrid(x1p, x2p)

Xp = np.column_stack((X1p.ravel(), X2p.ravel()))

yPred, yStd = gp.predict(Xp, return_std=True)

# Reshape predictions to match the grid shape (X1, X2)
yPred = np.asarray(yPred).reshape(X1p.shape)
yStd = np.asarray(yStd).reshape(X1p.shape)

plt.figure(2)
plt.contourf(X1p, X2p, yPred, cmap='jet', edgecolor='none')
plt.colorbar(label='y')
plt.xlabel('$X_1$')
plt.ylabel('$X_2$')  
plt.plot(X1p,X2p,'ko')
plt.subplots_adjust(bottom=0.18,  right=0.85,left=0.25)


plt.figure(3)
plt.contourf(X1p, X2p, yStd, cmap='jet', edgecolor='none')
plt.colorbar(label=r'$\sigma$')
plt.xlabel('$X_1$')
plt.ylabel('$X_2$')  
plt.subplots_adjust(bottom=0.18,  right=0.9,left=0.25)


plt.show()
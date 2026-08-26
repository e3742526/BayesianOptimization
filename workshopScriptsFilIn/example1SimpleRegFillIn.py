import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# Input data (padded with 1's so it = 1 + x)
X = np.array([[1, 1, 1],[-1,0,1]])
# output Data
y = np.array([1, 2, 4])

###################
# Fill In Code Here
# Prior (Signal) Covariance
s = 50
Sigma = np.array([[1, 0], [0, s]])
# Data noise  Variance
noise = 0.5
###################


# Interpoloation data points X*
xStar = np.zeros((2, 100))
xStar[0, :] = np.ones(100)
xStar[1, :] = np.linspace(X[1].min(), X[1].max(), 100)      

# check that the eigenvalue are positive (i.e. positive semi-definite)
print(np.linalg.eigvalsh(Sigma))

A = noise**-2 * (X @ X.T) + np.linalg.inv(Sigma)

# Mean prediction
mean = noise**-2 * xStar.T @ np.linalg.inv(A) @ X @ y

# Covariance prediction
stdevMat = xStar.T @ np.linalg.inv(A) @ xStar 
# standard devation from the covariane matrix
stdev = np.sqrt(np.diag(stdevMat))

matplotlib.rcParams.update({"font.size": 25})
plt.figure(1)
plt.plot(X[1,:].ravel(), y, 'ro')
plt.plot(xStar[1, :].ravel(), mean, 'b-')
plt.plot(xStar[1, :].ravel(), mean + stdev, 'b--')
plt.plot(xStar[1, :].ravel(), mean - stdev, 'b--')

plt.xlabel('X')
plt.ylabel('y')
plt.subplots_adjust(bottom=0.18, left=0.25)

plt.show()

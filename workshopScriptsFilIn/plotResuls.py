import numpy as np
import matplotlib.pyplot as plt

plotExperiment = True 

plt.rcParams.update({'font.size': 25})

data = np.loadtxt('tv.txt')
timeSim = data[:, 0]
velocitySim = 1000*data[:, 1]

dataData = np.loadtxt('TaData.txt',delimiter=',')
timeData = dataData[:, 0]
velocityData = dataData[:, 1]

plt.plot(timeSim, velocitySim)

if plotExperiment:  
    plt.plot(timeData, velocityData)

plt.xlabel('Time, $\\mu$s')
plt.ylabel('Velocity, m/s')
plt.subplots_adjust(left=0.20, bottom=0.20)

plt.show()
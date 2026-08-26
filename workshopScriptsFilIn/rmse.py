def rmse(simDataString, dataDataString):
    import numpy as np
    try:
        data = np.loadtxt(simDataString)
        dataData = np.loadtxt(dataDataString,delimiter=',')
    except:
        print('KO simulation data needs to be first argument and experimental data needs to be second argument')

    timeSim = data[:, 0]
    velocitySim = 1000*data[:, 1]
    startTimeSim = timeSim[0]
    endTimeSim = timeSim[-1]
    
    timeData = dataData[:, 0]
    velocityData = dataData[:, 1]
    startTimeData = timeData[0]
    endTimeData = timeData[-1]

    startTime = max(startTimeSim, startTimeData)
    endTime = min(endTimeSim, endTimeData)

    startIndexSim = np.searchsorted(timeSim, startTime, side='right')
    endIndexSim = np.searchsorted(timeSim, endTime, side='left')

    startIndexData = np.searchsorted(timeData, startTime, side='right')
    endIndexData = np.searchsorted(timeData, endTime, side='left')

    timeSim = timeSim[startIndexSim:endIndexSim]
    velocitySim = velocitySim[startIndexSim:endIndexSim]

    timeData = timeData[startIndexData:endIndexData]
    velocityData = velocityData[startIndexData:endIndexData]

    lengthSim = len(timeSim)
    lengthData = len(timeData)

    if lengthSim < lengthData:
        timeSmall = timeSim
        velocitySmall = velocitySim
        timeLarge = timeData
        velocityLarge = velocityData

    else:
        timeSmall = timeData
        velocitySmall = velocityData
        timeLarge = timeSim
        velocityLarge = velocitySim

    velocityLargeInterp = np.interp(timeSmall, timeLarge, velocityLarge)

    rmse = np.sqrt(np.mean((velocitySmall - velocityLargeInterp) ** 2))

    return rmse

if __name__ == "__main__":
    import sys

    simDataString = sys.argv[1] if len(sys.argv) > 1 else 'tv.txt'
    dataDataString = sys.argv[2] if len(sys.argv) > 2 else 'TaData.txt'

    result = rmse(simDataString, dataDataString)

    print('The RMSE is:', result)

import subprocess

cpuType = 'Mac'
#cpuType = 'PC'

# Compile KO with ko.in input deck paarameters
subprocess.run(['gfortran', 'KOv13.f'], check=True)

# run KO
if cpuType == 'PC':
    subprocess.run(['./a.exe'], check=True, text=True)
else:
    subprocess.run(['./a.out'], check=True, text=True)
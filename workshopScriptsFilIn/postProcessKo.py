import sys
import re
from io import StringIO
import numpy as np


def postProcessKo():
    with open("ko.dat") as f:
        text = f.read()

    # Insert an E before exponents that are missing one
    text = re.sub(
        r'([0-9]\.[0-9]+)([+-]\d{2,3})(?=\s|$)',
        r'\1E\2',
        text,
    )

    ko = np.loadtxt(StringIO(text))

    # =============================================================================
    # User settings
    # =============================================================================

    ix = 5      # MATLAB column 6 (0-based indexing)
    iy = 10     # MATLAB column 11

    tracer = 306 + 960 -1

    # Optional interface list (MATLAB node numbers)
    # interface = np.array([0, 594, 594+1710])

    # =============================================================================
    # Internal calculations
    # =============================================================================

    jj = int(ko[0, 2])

    if 'interface' in locals():
        interface = interface // 2

    itracer = tracer // 2 + 1
    itotal = jj // 2 - 1

    n = int(2 * len(ko) / (jj - 2))

    # =============================================================================
    # Determine scaling
    # =============================================================================

    if iy == 10:      # particle velocity
        convert = 10.0
        ystr = "Particle Velocity [km/s]"

    elif iy == 11:    # pressure
        convert = 100.0
        ystr = "Pressure [GPa]"

    elif iy == 19:    # temperature
        convert = 1.0
        ystr = "Temperature [K]"

    elif iy == 20:    # longitudinal stress
        convert = -100.0
        ystr = "Longitudinal Stress [GPa]"

    elif iy == 21:    # lateral stress
        convert = -100.0
        ystr = "Lateral Stress [GPa]"

    else:
        convert = 1.0
        ystr = ""


    velocity = np.zeros(n)
    time = np.zeros(n)

    for i in range(n):

        start = i * itotal
        end = (i + 1) * itotal

        row = start + itracer

        x = ko[row, ix]
        y = convert * ko[row, iy]

        velocity[i] = y
        time[i] = x

    # Write time and velocity to a text file
    tv = np.column_stack((time, velocity))

    np.savetxt(
        "tv.txt",
        tv,
        fmt="%.8e",
        comments=""
    )

if __name__ == "__main__":
    postProcessKo()


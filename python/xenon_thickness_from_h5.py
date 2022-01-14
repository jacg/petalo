"""Calculate thickness of LXe from data in HDF5 MC/vertices

Requires 1 command line argument: the name of the .h5 file
"""

import sys
import h5py
import numpy as np

filename = sys.argv[1]
vertices = h5py.File(filename)['MC']['vertices']
in_LXe = vertices.fields('volume_id')[:] == 0
x      = vertices.fields('x')[in_LXe]
y      = vertices.fields('y')[in_LXe]
r = np.sqrt(x*x + y*y)
r_max = r.max()
r_min = r.min()
LXe_thickness = r_max - r_min
print(f'{r_max:.1f} - {r_min:.1f} = {LXe_thickness:.1f} mm')

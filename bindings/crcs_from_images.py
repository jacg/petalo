import struct
from collections import namedtuple
from glob import glob
from math import sin, cos, pi
from fulano import FomConfig

def read_raw(filename):
    buffer = open(filename, 'rb').read()
    return struct.unpack_from('f' * (len(buffer) // 4), buffer)

def polar(r, phi):
    return (r * cos(phi), r * sin(phi))

# Helpers for defining the region of interest
cylinderX = namedtuple('cylinderX', 'y,z,r')
cylinderY = namedtuple('cylinderY', 'x,z,r')
cylinderZ = namedtuple('cylinderZ', 'x,y,r')
sphere    = namedtuple('sphere',  'x,y,z,r')
# Any python object with the thus-named attributes will do.

# Specify the locations of the regions of interest
step = pi / 6
roi_from_centre = 50
hot, cold, bg_radius, bg_activity = 4, 0, 4, 1

rois = ((sphere(*polar(roi_from_centre,  2*step), 0,  4.0),  hot),
        (sphere(*polar(roi_from_centre,  4*step), 0,  6.5),  hot),
        (sphere(*polar(roi_from_centre,  6*step), 0,  8.5),  hot),
        (sphere(*polar(roi_from_centre,  8*step), 0, 11.0),  hot),
        (sphere(*polar(roi_from_centre, 10*step), 0, 14.0), cold),
        (sphere(*polar(roi_from_centre, 12*step), 0, 18.5), cold))

bg_rois = (sphere(*polar(roi_from_centre,  1*step), 0, bg_radius),
           sphere(*polar(roi_from_centre,  3*step), 0, bg_radius),
           sphere(*polar(roi_from_centre,  5*step), 0, bg_radius),
           sphere(*polar(roi_from_centre,  7*step), 0, bg_radius),
           sphere(*polar(roi_from_centre,  9*step), 0, bg_radius),
           sphere(*polar(roi_from_centre, 11*step), 0, bg_radius))

# Specify the size and resolution of the field of view
size   = (180,) * 3
voxels =  (60,) * 3

# Encapsulate the above data
cfg = FomConfig(rois, bg_rois, bg_activity, voxels, size)

pattern = '../data/out/mlem/60_60_60_tof_100_*.raw'

for filename in sorted(glob(pattern)):
    data = read_raw(filename)
    for datum in cfg.crcs(data):
        print(f'{datum:8.2f}', end='')
    print()

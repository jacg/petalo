from sys   import argv
from utils import read_raw
import struct

_, nx,ny,nz, dx,dy,dz, *filenames = argv

nx,ny,nz = map(int  , (nx,ny,nz))
dx,dy,dz = map(float, (dx,dy,dz))

for name in filenames:
    in_ = open(name     , 'rb').read()
    out = open(name+'_h', 'wb')
    pixels = struct.pack('>HHH', nx, ny, nz)
    mm     = struct.pack('>fff', dx, dy, dz)
    out.write(pixels)
    out.write(mm)
    out.write(in_)

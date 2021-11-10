from sys   import argv
from utils import read_raw
import struct

def usage(argv):
    return f"""Add size information headers to 3D raw image files

    Usage:

       {argv[0]}  NX NY NZ  DX DY DZ  FILENAMES...

    NX, NY, NZ: number of voxels in each dimension
    DX, DY, DZ: full extent of FOV in each dimension

    The original files will not be modified. The headers will be added to new
    files matching the original names with the suffix `_h`"""

if '-h' in argv or '--help' in argv:
    exit(usage(argv))

try:
    prog, nx,ny,nz, dx,dy,dz, *filenames = argv
except ValueError:
    exit(usage(argv))

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

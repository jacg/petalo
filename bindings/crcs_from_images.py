import struct
from glob import glob
from fulano import crcs

def read_raw(filename):
    buffer = open(filename, 'rb').read()
    return struct.unpack_from('f' * (len(buffer) // 4), buffer)

pattern = '../data/out/mlem/60_60_60_tof_100_*.raw'

for filename in sorted(glob(pattern)):
    data = read_raw(filename)
    for datum in crcs(data):
        print(f'{datum:8.2f}', end='')
    print()

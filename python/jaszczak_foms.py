#!/usr/bin/env python3

"""Calculate FOMs for Jaszczak image files in a directory

Usage: jaszczak_foms.py DIR

Arguments:
  DIR directory containing *NN.raw image files of Jaszczak phantom
"""

import sys
import re
from pathlib import Path
import subprocess
from collections import namedtuple
from matplotlib import pyplot as plt
from docopt import docopt


FOMS = namedtuple('foms', 'crcs, snrs')
CRC  = namedtuple('crc', 'crc, var')

xxx = dict(true='Truelor', reco='mlor')
TRUE, RECO = 'Truemlor', 'mlor'

sphere_diameters = 9.5, 12.7, 15.9, 19.1, 25.4, 31.8

def get_foms(image_file_name):

    command = f'cargo run --release --bin foms -- {image_file_name} jaszczak'
    o = subprocess.run(command, shell=True, capture_output=True)
    crcs = {}
    snrs = {}

    for line in o.stdout.decode('utf-8').split('\n'):

        columns = line.split()
        if len(columns) == 4:
            d, crc, var, snr = map(float, columns)
            crcs[d] = CRC(crc=crc, var=var)
            snrs[d] = snr

    return FOMS(crcs=crcs, snrs=snrs)


from contextlib import contextmanager
@contextmanager
def dummy(*args, **kwds):
    yield

cli_args = docopt(__doc__)

def write(*args, **kwds):
    print(*args, **kwds, file=sys.stdout)
    print(*args, **kwds, file=outfile)

directory = cli_args['DIR']
print(f'Calculating FOMs for images in {directory}')
filename = f'{directory}/foms'
images = sorted(Path(directory).glob('*.raw'))
with open(filename, 'w') as outfile:
    write('  ', end='')
    ds_header = ''.join(f'{d:7.1f}' for d in sphere_diameters)
    write(f'{"CRCs":>25}        {"background variabilities":>53}          {"SNRs":>31}', end='\n\n')
    write(f'  {ds_header}         {ds_header}         {ds_header}')
    write()
    for image_file_name in images:
        iteration = int(re.search(r'.*([0-9]{2}).raw$', str(image_file_name)).group(1))
        foms = get_foms(image_file_name)
        crcs, snrs = foms
        write(f'{iteration+1:2} ', end='')
        crcs, variabilities = tuple(zip(*((crc, var) for (r, (crc, var)) in crcs.items())))
        write(''.join(f'{c:6.1f} ' for c in crcs)         , end='         ')
        write(''.join(f'{v:6.1f} ' for v in variabilities), end='         ')
        write(''.join(f'{r:6.1f} ' for r in snrs.values())) # look broken in the Rust implementation of foms


#TODO add ability to read foms from generated fom files

exit(0)

fs = {}

run_numbers_to_use = (20,)

for run_number in run_numbers_to_use:
    for it in iterations:
        run_spec = runs[run_number]
        img_spec = run_spec + (it,)
        print(img_spec)
        print(get_foms(*img_spec))
        fs[img_spec] = get_foms(*img_spec)

for d in sphere_diameters:
    plt.figure()
    count = 0
    for run_number in run_numbers_to_use:
        run_spec = runs[run_number]
        img_spec = run_spec + (it,)
        y = [fs[img_spec].crcs[d].crc for i in iterations]
        e = [fs[img_spec].crcs[d].var for i in iterations]
        x = [1 + n - count/10 for n in iterations]
        #x = iterations
        #plt.errorbar(x, y, yerr=e, label=f'TOF={t} {c}', capsize=3)
        plt.errorbar(x, y, yerr=e, label=f'TOF=t c', capsize=3)
        count += 1
    plt.legend()
    plt.ylim(top=120)
    plt.title(f'CRC vs iteration for {d}mm sphere')
    plt.savefig(f'crcs-{d}mm.png')

plt.show()

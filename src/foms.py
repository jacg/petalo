#!/usr/bin/env python3

"""Calculate FOMs for jaszczak and nema7 image files in a directory

Usage:
    foms.py <phantom> <DIR> [options]

Arguments:
    phantom     #choose phantom: nema7 or jaszczak
    DIR         #directory containing *.raw image files of the specified phantom

options:
    --show-plot=<bool>      Display the plot generated for the FOMs [default=False]
    --hide-output=<bool>    Don't print results to stdout [default=False]
"""

import sys
import re
from pathlib import Path
import subprocess
from collections import namedtuple, defaultdict
from matplotlib import pyplot as plt
from docopt import docopt
from itertools import islice

CRC  = namedtuple('crc', 'crc, var')
NFOMS = namedtuple('foms', 'crcs, aocs')
JFOMS = namedtuple('foms', 'crcs, snrs')

sphere_diameters = dict(
    jaszczak = ( 9.5, 12.7, 15.9, 19.1, 25.4, 31.8),
    nema7    = (10.0, 13.0, 17.0, 22.0, 28.0, 37.0),
)

def get_foms(command):

    o = subprocess.run(command, shell=True, capture_output=True)

    crcs = {}
    snrs = {}

    for line in o.stdout.decode('utf-8').split('\n'):

        columns = line.split()
        if len(columns) == 4:
            d, crc, var, snr = map(float, columns)
            crcs[d] = CRC(crc=crc, var=var)
            snrs[d] = snr

        if line.startswith('AOCs: [') and line.endswith(']'):
            nocommas = ''.join(c for c in line[7:-1] if c != ',')
            aocs = tuple(map(float, nocommas.split()))

    return JFOMS(crcs=crcs, snrs=snrs)

def write(outfile, cli_args, *args, **kwds):
    if bool(cli_args['--hide-output']): 
        print(*args, **kwds, file=outfile)
    else:
        print(*args, **kwds, file=outfile)
        print(*args, **kwds, file=sys.stdout)

def write_command(filename, images, commandStart, sphere_diameters,cli_args):
    if Path(filename).is_file():
        with open(filename, 'r') as infile:
            for line in infile:
                print (line.rstrip())
        return
    with open(filename, 'w') as outfile:
        write(outfile,cli_args, '  ', end='')
        ds_header = ''.join(f'{d:7.1f}' for d in sphere_diameters)
        write(outfile,cli_args, f'{"CRCs":>25}        {"background variabilities":>53}          {"SNRs":>31}', end='\n\n')
        write(outfile,cli_args, f'     {ds_header}         {ds_header}         {ds_header}')
        write(outfile,cli_args, )
        for image_file_name in images:
            matches = re.search(r'.*([0-9]{2})-([0-9]{2}).raw$', str(image_file_name))
            iteration = int(matches.group(1))
            subset    = int(matches.group(2))
            foms = get_foms(f'{commandStart} {image_file_name}')
            crcs, snrs = foms
            write(outfile,cli_args, f'{iteration:02}-{subset:02} ', end='')
            crcs, variabilities = tuple(zip(*((crc, var) for (r, (crc, var)) in crcs.items())))
            write(outfile,cli_args, ''.join(f'{c:6.1f} ' for c in crcs)         , end='         ')
            write(outfile,cli_args, ''.join(f'{v:6.1f} ' for v in variabilities), end='         ')
            write(outfile,cli_args, ''.join(f'{r:6.1f} ' for r in snrs.values())) # look broken in the Rust implementation of foms
def plot_from_fom(dir, sphere_diameters, cli_args):

    iterations = []
    subsets    = []
    crcs = defaultdict(list)
    bgvs = defaultdict(list)
    snrs = defaultdict(list)

    with open(dir, encoding='utf-8') as f:
        for line in islice(f, 4, None):

            iteration_subset, *values = line.split()
            iteration, subset = map(int, iteration_subset.split('-'))
            crcs_ = map(float, values[  : 6])
            bgvs_ = map(float, values[ 6:12])
            snrs_ = map(float, values[12:  ])

            iterations.append(iteration)
            subsets   .append(subset)

            for d, crc, bgv, snr in zip(sphere_diameters, crcs_, bgvs_, snrs_):
                crcs[d].append(crc)
                bgvs[d].append(bgv)
                snrs[d].append(snr)

    for d in sphere_diameters:
        y = crcs[d]
        e = bgvs[d]
        x = tuple(range(len(y)))
        plt.figure()
        plt.plot(x,y,linewidth=2.0)
        plt.errorbar(x,y,yerr=e,label=f'TOF=t c',capsize=3)
        plt.legend()
        plt.ylim(top=120)
        plt.title(f'CRC vs iteration for {d}mm sphere')
        plt.savefig(f'{cli_args["<DIR>"]}/crcs-{d}mm.png')

    if cli_args['--show-plot']:
        plt.show()

def main():
    cli_args = docopt(__doc__)
    directory = cli_args['<DIR>']
    print(f'Calculating FOMs for images in {directory}')
    filename = f'{directory}/foms'
    images = sorted(Path(directory).glob('*.raw'))

    known_phantoms = ('jaszczak', 'nema7')
    phantom = cli_args['<phantom>']
    if phantom not in known_phantoms:
        sys.exit(f"You must specify a known phantom type: (jaszczak/nema7), not {phantom}")

    command = f'target/release/foms {phantom}'
    write_command(filename, images, command, sphere_diameters[phantom], cli_args)
    plot_from_fom(filename,                  sphere_diameters[phantom], cli_args)

main()

exit(0)

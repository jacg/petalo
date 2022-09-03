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
import numpy as np

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

def write_command(directory, images, commandStart, sphere_diameters,cli_args):
    filename = f'{directory}/foms'
    if Path(filename).is_file():
        with open(filename, 'r') as infile:
            for line in infile:
                print (line.rstrip())
        print("\nThese figures were retrieved from the cache file.")
        print(f"To recalculate the figures delete the cache file: {filename}")
        return

    subprocess.run("cargo build --release --bin foms", shell=True, capture_output=False)

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

def plot_from_fom(directory, sphere_diameters, cli_args):

    cache_file = f'{directory}/foms'

    iterations = []
    subsets    = []
    crcs = defaultdict(list)
    bgvs = defaultdict(list)
    snrs = defaultdict(list)

    with open(cache_file, encoding='utf-8') as f:
        for line in islice(f, 4, None):

            iteration_subset, *values = line.split()
            iteration, subset = map(int, iteration_subset.split('-'))
            crcs_ = map(float, values[  : 6])
            bgvs_ = map(float, values[ 6:12])
            snrs_ = map(float, values[12:  ])

            iterations.append(iteration)
            subsets   .append(subset)

            for d, crc, bgv, snr in zip(sphere_diameters, crcs_, bgvs_, snrs_):
                crcs[d].append(crc / 100)
                bgvs[d].append(bgv / 100)
                snrs[d].append(snr / 100)

    plt.rc('font'  ,     size=20)
    plt.rc('legend', fontsize=10)

    fig, ax_crc = plt.subplots(figsize=(16, 12))
    ax_crc.plot()
    ax_snr = ax_crc.twinx()
    ax_crc.grid(axis='y')
    ax_crc.set_yticks(np.linspace(0, 1,11))
    ax_snr.set_yticks(np.linspace(0,10,11))

    # Show which linestyle corresponds to which FOM, in legend.
    c_style = '-'
    s_style = '--'
    b_style = (0, (10,3))
    ax_crc.plot([], [], color='black', linestyle=c_style, label='CRC')
    ax_crc.plot([], [], color='black', linestyle=s_style, label='SNR')
    ax_crc.plot([], [], color='black', linestyle=b_style, label='bg-var')

    n_subsets    = max(subsets)
    n_iterations = max(iterations)
    colors = 'blue orange purple green cyan red'.split()
    for d, color in reversed(tuple(zip(sphere_diameters, colors))):
        c = crcs[d]
        s = snrs[d]
        b = bgvs[d]
        x = tuple(range(len(c)))

        # Label axes differently for OSEM vs MLEM
        if all(subset == 1 for subset in subsets):
            # MLEM: Plain iteration numbers
            x = tuple(range(len(c)))
            ax_crc.set_xlabel('iteration')
        else:
            # OSEM: iteration_number:subset_number
            x = [f'{i}:{s}' for i,s in zip(iterations, subsets)]
            ax_crc.set_xlabel('iteration : subset')

        ax_crc.plot(x, c, color=color, linestyle=c_style, marker=' ', linewidth=2.0, label=f'{d}mm')
        ax_snr.plot(x, s, color=color, linestyle=s_style, marker=' ', linewidth=2.0, label=None)
        ax_crc.plot(x, b, color=color, linestyle=b_style, marker=' ', linewidth=2.0, label=None)
        #plt.errorbar(x,y,yerr=e,label=f'{d}mm',capsize=3)
    ax_crc.set_ylim(bottom=0, top= 1); ax_crc.set_ylabel('CRC / bg-var')
    ax_snr.set_ylim(bottom=0, top=10); ax_snr.set_ylabel('SNR')
    ax_crc.legend()

    # Vertical lines separating OSEM iterations
    if n_subsets > 1:
        # OSEM (not MLEM)
        for n, (i, s) in enumerate(zip(iterations, subsets)):
            if n == 0: continue
            if s == 1:
                ax_crc.plot([n, n], [0,1], color='gray', linewidth=0.5)
        ax_crc.set_xticks(tuple(n-1 for n in range(0, len(x)+1, 5)))
        title = 'CRCs and SNRs vs iteration:subset'
    else:
        # MLEM (not OSEM)
        title = 'CRCs and SNRs vs iteration'
    plt.title(title)

    suptitle = get_plot_title(directory)
    if suptitle:
        plt.suptitle(suptitle.strip())

    plt.savefig(f'{cli_args["<DIR>"]}/foms.png', dpi=80)

    if cli_args['--show-plot']:
        plt.show()


def get_plot_title(directory):
    path = Path(f'{directory}/title')
    if path.is_file():
        return path.read_text()


def main():
    cli_args = docopt(__doc__)
    directory = cli_args['<DIR>']
    print(f'Calculating FOMs for images in {directory}')
    images = sorted(Path(directory).glob('*.raw'))

    known_phantoms = ('jaszczak', 'nema7')
    phantom = cli_args['<phantom>']
    if phantom not in known_phantoms:
        sys.exit(f"You must specify a known phantom type: (jaszczak/nema7), not {phantom}")

    command = f'target/release/foms {phantom}'
    write_command(directory, images, command, sphere_diameters[phantom], cli_args)
    plot_from_fom(directory,                  sphere_diameters[phantom], cli_args)


main()

#!/usr/bin/env python3

"""Calculate FOMs for jaszczak and nema7 image files in a directory

Usage:
    foms.py <jaszczak> <DIR> [options]
    foms.py <nema7> <DIR> [options]

Arguments:
    jaszczak    #specifies to work with a jaszczak type phantom
    nema7       #specifies to work with a nema7 type phantom
    DIR         #directory containing *.raw image files of the specified phantom

options:
    --show-plot=<bool>      Display the plot generated for the FOMs [default=False]
    --hide-output=<bool>    Don't print results to stdout [default=False]
"""

import sys
import re
from pathlib import Path
import subprocess
from collections import namedtuple
from matplotlib import pyplot as plt
from docopt import docopt

CRC  = namedtuple('crc', 'crc, var')
NFOMS = namedtuple('foms', 'crcs, aocs')
JFOMS = namedtuple('foms', 'crcs, snrs')

nema7_sphere_diameters = 10.0, 13.0, 17.0, 22.0, 28.0, 37.0
jaszczak_sphere_diameters = 9.5, 12.7, 15.9, 19.1, 25.4, 31.8

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
            foms = get_foms(commandStart + str(image_file_name))
            crcs, snrs = foms
            write(outfile,cli_args, f'{iteration:02}-{subset:02} ', end='')
            crcs, variabilities = tuple(zip(*((crc, var) for (r, (crc, var)) in crcs.items())))
            write(outfile,cli_args, ''.join(f'{c:6.1f} ' for c in crcs)         , end='         ')
            write(outfile,cli_args, ''.join(f'{v:6.1f} ' for v in variabilities), end='         ')
            write(outfile,cli_args, ''.join(f'{r:6.1f} ' for r in snrs.values())) # look broken in the Rust implementation of foms
    
def plot_from_fom(dir, sphere_diameters, cli_args):
    count = 0
    isData = False
    data = dict()
    for d in sphere_diameters:
        data[d] = [] #initialize lists for each diameter key
    with open(dir, encoding='utf-8') as f:
        Lines = f.readlines()
        for line in Lines:
            if count == 4:
                isData = True
            if not isData:
                count += 1
                continue
            #remove first element from a line and convert the remaining elements to floats
            parsedLine = list(map(float, line.split()[1:]))
            for i in range(len(sphere_diameters)):
                data[sphere_diameters[i]].append((float(count-4), parsedLine[i],
                                                  parsedLine[i+len(sphere_diameters)]))
            count += 1
    for d in sphere_diameters:
        x = [data[d][i][0] for i in range(len(data[d]))]
        y = [data[d][i][1] for i in range(len(data[d]))]
        e = [data[d][i][2] for i in range(len(data[d]))]
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

    if bool(cli_args['<jaszczak>']):
        commandStart = 'target/release/foms -- jaszczak '
        write_command(filename, images, commandStart, jaszczak_sphere_diameters,
                      cli_args)
        plot_from_fom(filename, jaszczak_sphere_diameters, cli_args)

    elif bool(cli_args['<nema7>']):
        commandStart = 'target/release/foms -- nema7 '
        write_command(filename, images, commandStart, nema7_sphere_diameters,
                      cli_args)
        plot_from_fom(filename, nema7_sphere_diameters, cli_args)

    else:
        print("You must specify a known phantom type: (jaszczak/nema7)")

if __name__=='__main__':
   main()

#TODO add ability to read foms from generated fom files

            
    

"""
def plot(sphere_diameters, iterations, runs):  
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
"""

exit(0)

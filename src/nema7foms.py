#!/usr/bin/env python3

import subprocess
from collections import namedtuple
from matplotlib import pyplot as plt

FOMS = namedtuple('foms', 'crcs, aocs')
CRC  = namedtuple('crc', 'crc, var')

tofs = 'NO', '250'
cors = 'uncorrected', # 'corrected' #, 'corrected-flipped'
iterations = (0, 4, 9, 14, 19, 24, 29)
iterations = tuple(range(30))

file_base = 'data/nema7/images/1m-steel-signal-lors-true-again'
file_base = 'data/nema7/images/nema7-airbody-dz1m-LXe40mm-noQrz-5-lors-true'

def get_foms(tof, corr, it):

    command = f'cargo run --release --bin nema7_analysis -- {file_base}-tof{tof}-{corr}_{it:02}.raw'
    o = subprocess.run(command, shell=True, capture_output=True)
    crcs = {}

    for line in o.stdout.decode('utf-8').split('\n'):


        columns = line.split()
        if len(columns) == 3:
            crc, var, d = columns
            d = d[1:-1]
            crc = float(crc)
            var = float(var)
            d   =   int(d)
            crcs[d] = CRC(crc=crc, var=var)

        if line.startswith('[') and line.endswith(']'):
            nocommas = ''.join(c for c in line[1:-1] if c != ',')
            aocs = tuple(map(float, nocommas.split()))

    return FOMS(crcs=crcs, aocs=aocs)


fs = {}

for it in iterations:
    for tof in tofs:
        for corr in cors:
            fs[(tof, corr, it)] = get_foms(tof, corr, it)

for d in (10, 13, 17, 22, 28, 37):
    plt.figure()
    count = 0
    for t in tofs:
        for c in cors:
            y = [fs[t,c,i].crcs[d].crc for i in iterations]
            e = [fs[t,c,i].crcs[d].var for i in iterations]
            x = [1 + n - count/10 for n in iterations]
            #x = iterations
            plt.errorbar(x, y, yerr=e, label=f'TOF={t} {c}', capsize=3)
            count += 1
    plt.legend()
    plt.ylim(top=120)
    plt.title(f'CRC vs iteration for {d}mm sphere')
    plt.savefig(f'crcs-{d}mm.png')

plt.show()

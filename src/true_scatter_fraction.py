#!/usr/bin/env python3

import argparse
from operator import itemgetter
from itertools import groupby, chain
import h5py

parser = argparse.ArgumentParser()
parser.add_argument('infiles', nargs='+', help='HDF5 input files with MC events')
args = parser.parse_args()

lxe_vol = 0

n_files = len(args.infiles)

def vertices_in_file(filename):
    return iter(h5py.File(filename)['MC']['vertices'])


def in_lxe_or_cylinder(vertex):
    return vertex[-1] in (lxe_vol, cylinder_vol)

def in_lxe (vertex): return vertex[-1] == lxe_vol
def primary(vertex): return vertex[ 1] <= 2
def trackid(vertex): return vertex[ 1]
def volume (vertex): return vertex[-1]
def E      (vertex): return vertex[ 8]

lor_like = 0
scatter  = 0

for i, filename in enumerate(args.infiles, 1):
    vertices = vertices_in_file(filename)
    vertices_in_xe = filter(in_lxe, vertices)
    grouped_vertices = groupby(vertices_in_xe, itemgetter(0))

    for eventid, vertices_in_evt in grouped_vertices:
        first_enegies = {tid:E(next(vertices))
                         for tid, vertices in groupby(vertices_in_evt, trackid)
                         if tid <= 2}
        if len(first_enegies) == 2:
            lor_like += 1
            if min(first_enegies.values()) < 511:
                scatter += 1
    print(f'{scatter} / {lor_like} = {100 * scatter / lor_like:.1f}% ' +
          f' scatter fraction after {i:4} / {n_files} files')

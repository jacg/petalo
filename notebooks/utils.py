import struct
from collections import namedtuple
import os.path
import subprocess

import numpy as np
from operator import mul
from functools import reduce

import matplotlib.pyplot as plt


#def compare_c_rust(*, tof=None, shape=(60,60,60), iterations=6, use_true=False, file_in=None, n_events):
def compare_c_rust(**incoming):
    # default values
    args = dict(tof=None, shape=(60,60,60), size=(360, 360, 360),
                iterations=6, use_true=False, file_in=None, c=False,
                dataset='', n_events=False)
    # override defaults with incoming values
    args.update(incoming)
    #args = Args(tof, shape, iterations, c=False, use_true=use_true, file_in=file_in, n_events=n_events)
    args = Args(**args)
    print(args)
    display_reconstructed_images(nr=1, nc=6, args=args)
    display_reconstructed_images(nr=1, nc=6, args=update_nt(args, c=True))


def display_reconstructed_images(nr, nc, args):
    generate_data_if_missing(args)
    fig, ax = plt.subplots(nr,nc, figsize=(nc*5, nr*5))
    fig.suptitle("C" if args.c else "Rust")

    for (r,c) in ((r,c) for r in range(nr) for c in range(nc)):
        filename = args_to_filename(args, c+nc*r)
        the_ax = ax[c] if min(nr,nc) == 1 else ax[r,c]
        show_image_from_file(filename, args, ax=the_ax)


def show_image_from_file(filename, args, ax=plt):
    data = read_raw(filename)
    image = wrap_1d_into_3d(data, args.shape)
    show_image(image, zslice=args.shape[-1]//2, extent=args.size, ax=ax)


def show_image(image, zslice, extent, ax=plt):
    xe, ye, ze = extent
    ax.imshow(image[:, :, zslice].transpose(),
              extent = [-xe,xe, -ye,ye],
              origin = 'lower')


def read_raw(filename):
    buffer = open(filename, 'rb').read()
    return struct.unpack_from('f' * (len(buffer) // 4), buffer)


def wrap_1d_into_3d(data, shape, row_major=False):
    assert len(data) == reduce(mul, shape)
    image = np.zeros(shape)
    nx, ny, nz = shape
    for index, x in enumerate(data):
        i = int(index % nx)
        j = int(index / nx) % ny
        k = int(index / (nx * ny))
        if row_major: image[k,j,i] = x
        else        : image[i,j,k] = x
    return image


def update_nt(nt, **overrides):
    return type(nt)(**{**nt._asdict(), **overrides})


def generate_data_if_missing(args):
    required = args_to_filename(args)
    if os.path.exists(required):
        return
    command = args_to_cli(args)
    print("The required image files are missing. Generating them with:\n")
    print(f'   {command}\n')
    original_dir = os.getcwd()
    os.chdir('..')
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print('============================== FAILED: stderr ==================================')
        print(result.stderr)
        print('================================================================================')
    os.chdir(original_dir)


Args = namedtuple('Args', '''tof, shape, iterations, c, use_true, file_in, n_events, dataset,
                             read_lors, size''')


def args_to_cli(args):
    nx, ny, nz = args.shape
    sx, sy, sz = args.size
    n = f' -n {nx},{ny},{nz}'     if args.shape != (60,60,60) else ''
    s = f' -s {sx},{sy},{sz}'     if args.size                else ''
    i = f' -i {args.iterations}'
    c =  ' -c'                    if args.c                   else ''
    r = f' -r {args.tof}'         if args.tof                 else ''
    t =  ' --use-true'            if args.use_true            else ''
    f = f' -f {args.file_in}'     if args.file_in             else ''
    e = f' -e 0..{args.n_events}' if args.n_events            else ''
    d = f' -d {args.dataset}'     if args.dataset             else ''
    l = f' --read-lors'           if args.read_lors           else ''
    return f'cargo run --bin mlem --release -- {i}{r}{n}{s}{c}{t}{f}{e}{d}{l}'


def args_to_filename(args, N=None):
    nx, ny, nz = args.shape
    tof = args.tof
    TOF = tof if tof else 'OFF'
    c = 'c' if args.c else ''
    n = N if N is not None else (args.iterations - 1)
    return f'../data/out/{c}mlem/{nx}_{ny}_{nz}_tof_{TOF}_{n:02}.raw'

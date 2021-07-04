import struct
from collections import namedtuple
import os.path
import subprocess

import numpy as np
from operator import mul
from functools import reduce

import matplotlib.pyplot as plt


def compare_c_rust(**incoming):
    # default values
    args = dict(tof=None, shape=(60,60,60), full_length=(180, 180, 180),
                iterations=6, use_true=False, file_in=None, c=False,
                dataset='', n_events=False, axis = 'z', slice_='half')
    # override defaults with incoming values
    args.update(incoming)
    args = Args(**args)
    print(args)
    # Always run the Rust version, only run C if requested
    display_reconstructed_images(nr=1, nc=6, args=update_nt(args, c=False))
    if args.c:
        display_reconstructed_images(nr=1, nc=6, args=args)


def display_reconstructed_images(nr, nc, args):
    generate_data_if_missing(args)
    fig, ax = plt.subplots(nr,nc, figsize=(nc*5, nr*5))
    fig.suptitle("C" if args.c else "Rust")

    for (r,c) in ((r,c) for r in range(nr) for c in range(nc)):
        filename = args_to_filename(args, c+nc*r)
        the_ax = ax[c] if min(nr,nc) == 1 else ax[r,c]
        show_image_from_file(filename, axis=args.axis, slice_=args.slice_, ax=the_ax)


def load_image_from_file(filename):
    (shape, full_extent), data = read_raw(filename)
    image = wrap_1d_into_3d(data, shape)
    return image, full_extent


def show_image_from_file(filename, *, axis, slice_, ax=plt):
    (shape, full_length), data = read_raw(filename)
    image = wrap_1d_into_3d(data, shape)
    show_image(image, axis=axis, slice_=slice_, extent=full_length, shape=shape, ax=ax)


def show_image(image, *, axis, slice_, extent, shape, ax=plt):
    xe, ye, ze = (e/2 for e in extent)

    if axis == 'x': e, n = xe, shape[0]
    if axis == 'y': e, n = ye, shape[1]
    if axis == 'z': e, n = ze, shape[2]

    s = int(n * (slice_ / (2*e) - 0.5))

    if axis == 'x': ax.imshow(image[s, :, :]   , extent = [-ze,ze, -ye,ye], origin = 'lower')
    if axis == 'y': ax.imshow(image[:, s, :]   , extent = [-ze,ze, -xe,xe], origin = 'lower')
    if axis == 'z': ax.imshow(image[:, :, s].T , extent = [-xe,xe, -ye,ye], origin = 'lower')


def read_raw(filename):
    buffer = open(filename, 'rb').read()
    metadata_length = 18
    metadata = buffer[: metadata_length  ]
    data     = buffer[  metadata_length :]
    pixels = struct.unpack_from('>HHH', metadata[:6])
    mm     = struct.unpack_from('>fff', metadata[6:])
    data   = struct.unpack_from('>'+'f' * (len(data) // 4), data)
    return (pixels, mm), data


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
                             read_lors, full_length, axis, slice_''')


def args_to_cli(args):
    nx, ny, nz = args.shape
    sx, sy, sz = args.full_length
    n = f' -n {nx},{ny},{nz}'     if args.shape != (60,60,60) else ''
    s = f' -s {sx},{sy},{sz}'     if args.full_length         else ''
    i = f' -i {args.iterations}'
    c =  ' -c'                    if args.c                   else ''
    uc = ' --features ccmlem'     if args.c                   else ''
    r = f' -r {args.tof}'         if args.tof                 else ''
    t =  ' --use-true'            if args.use_true            else ''
    f = f' -f {args.file_in}'     if args.file_in             else ''
    e = f' -e 0..{args.n_events}' if args.n_events            else ''
    d = f' -d {args.dataset}'     if args.dataset             else ''
    l = f' --read-lors'           if args.read_lors           else ''
    return f'cargo run --bin mlem --release {uc} -- {i}{r}{n}{s}{c}{t}{f}{e}{d}{l}'


def args_to_filename(args, N=None):
    nx, ny, nz = args.shape
    tof = args.tof
    TOF = tof if tof else 'OFF'
    c = 'c' if args.c else ''
    n = N if N is not None else (args.iterations - 1)
    return f'../data/out/{c}mlem/{nx}_{ny}_{nz}_tof_{TOF}_{n:02}.raw'

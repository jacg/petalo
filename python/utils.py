import struct
from collections import namedtuple
import os.path
import subprocess
from operator import itemgetter

import numpy as np
from operator import mul
from functools import reduce

import matplotlib.pyplot as plt


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


def read_raw(filename, header_only=False):
    buffer = open(filename, 'rb').read()
    metadata_length = 18
    metadata = buffer[: metadata_length  ]
    data     = buffer[  metadata_length :]
    pixels = struct.unpack_from('>HHH', metadata[:6])
    mm     = struct.unpack_from('>fff', metadata[6:])
    if header_only:
        return pixels, mm
    data   = struct.unpack_from('>'+'f' * (len(data) // 4), data)
    return (pixels, mm), data


def wrap_1d_into_3d(data, shape, row_major=False):
    size_data = len(data)
    size_expected = reduce(mul, shape)
    if size_data != size_expected:
        x,y,z = shape
        exit(f'Error: Data length {size_data} does not match voxelization {x} x {y} x {z} = {size_expected}')
    image = np.zeros(shape)
    nx, ny, nz = shape
    for index, x in enumerate(data):
        i = int(index % nx)
        j = int(index / nx) % ny
        k = int(index / (nx * ny))
        if row_major: image[k,j,i] = x
        else        : image[i,j,k] = x
    return image


def read_raw_without_header(filename):
    data = open(filename, 'rb').read()
    data = struct.unpack_from('>'+'f' * (len(data) // 4), data)
    return data

from sys import argv
import argparse
from itertools import starmap

import numpy as np
import matplotlib.pyplot as plt

from sklearn.cluster import KMeans

import utils

parser = argparse.ArgumentParser()
parser.add_argument('infile', type=str  , help="Raw image file (containing size metadata) to be analysed")
parser.add_argument('--peak-box', type=float, help='box width around peak', default=25)
args = parser.parse_args()

def voxel_centres(image, full_extent):
    def coordinates_of_voxel_centres(n, dx):
        last_centre = (n-1) / n * dx / 2
        return np.linspace(-last_centre, last_centre, num=n)
    x,y,z = starmap(coordinates_of_voxel_centres, zip(image.shape, full_extent))
    coords = np.swapaxes(np.array(np.meshgrid(z,y,x)), 0,3)
    return np.flip(coords, axis=3)

image, full_extent = utils.load_image_from_file(args.infile)

positions_of_active_voxels = voxel_centres(image, full_extent)[image > 1]

kmeans = KMeans(n_clusters=6, random_state=0).fit(positions_of_active_voxels)

def zone(x):
    if x > 150: return 2
    if x >  50: return 1
    return 0

def by_zone(xyz): return tuple(map(zone, xyz))

cluster_centres = sorted(kmeans.cluster_centers_, key=by_zone)

def slice_to_position(slice_number, n_pixels, full_extent):
    pixel_width = full_extent / n_pixels
    return pixel_width * (slice_number + 0.5 - n_pixels / 2)

def position_to_slice(position, n_pixels, full_extent):
    return round((n_pixels * position) / full_extent + 0.5 * (n_pixels - 1))

def test_roundtrip(slice_number, n_pixels, full_extent):
    position  = slice_to_position(slice_number, n_pixels, full_extent)
    new_slice = position_to_slice(position    , n_pixels, full_extent)
    assert new_slice == slice_number


delta = args.peak_box
ex, ey, ez = full_extent
nx, ny, nz = image.shape
dx, dy, dz = (e/n for e,n in zip(full_extent, image.shape))

def analyse_one_point(x,y,z):

    ylo = position_to_slice(y-delta, ny, ey)
    yhi = position_to_slice(y+delta, ny, ey) + 1

    zlo = position_to_slice(z-delta, nz, ez)
    zhi = position_to_slice(z+delta, nz, ez) + 1

    around_peak = image[:, ylo:yhi, zlo:zhi]
    xpk, ypk, zpk = np.unravel_index(np.argmax(around_peak), around_peak.shape)
    peak_value = around_peak[xpk, ypk, zpk]

    fig, axs = plt.subplots(1,3, tight_layout=True)

    along_x = around_peak[ : , ypk, zpk]
    along_y = around_peak[xpk,  : , zpk]
    along_z = around_peak[xpk, ypk,  : ]

    print(f'{x:6.1f} {y:6.1f} {z:6.1f}     ', end='')

    widths = []
    for i, (a, d) in enumerate(zip((along_x, along_y, along_z), (dx, dy, dz))):
        hm = peak_value /  2
        tm = peak_value / 10
        lh, rh = full_width_at(hm, a, d)
        lt, rt = full_width_at(tm, a, d)
        ax = axs[i]
        ax.plot(a, '-o')
        ax.plot((lh, rh), (hm, hm))
        ax.plot((lt, rt), (tm, tm))
        fig.suptitle(f'{x:6.1f} {y:6.1f} {z:6.1f}')
        ax.set_title(f'FWHM = {rh-lh:4.1f} mm   FWTM = {rt-lt:4.1f} mm')
        widths.append((rh-lh, rt-lt))

    (hx,tx), (hy,ty), (hz,tz) = widths
    print(f'   {hx:4.1f} {hy:4.1f} {hz:4.1f}     {tx:4.1f} {ty:4.1f} {tz:4.1f}')


def n_pixels_to(height, array):
    for i, (this, prev) in enumerate(zip(array[1:], array), 1):
        if this > height:
            return i - 0.5 + (height - prev) / (this - prev)

def full_width_at(height, array, pixel_size):
    l = n_pixels_to(height, array)
    r = n_pixels_to(height, array[::-1])
    return (l - 0.5) * pixel_size, (len(array) - r - 0.5) * pixel_size

print('       centroid                FWHM / mm          FWTM / mm')
print('                              x    y    z        x    y    z')
for c in cluster_centres:
    analyse_one_point(*c)
plt.show()

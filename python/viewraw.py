"""View 3D raw image files

Usage: viewraw.py [--little-endian]                                     FILES...
       viewraw.py --show-header                                         FILES...
       viewraw.py [--assume-header NX NY NZ DX DY DZ] [--little-endian] FILES...
       viewraw.py    [--add-header NX NY NZ DX DY DZ]                   FILES...

Arguments:
  FILES  3D raw image files to be loaded into viewer

Options:
  --assume-header  View headerless raw 3D image, specifying number of voxels
                   and physicals size of FOV
  --little-endian  Assume the input files contain little-endian data (rather
                   than big-endian, as guaranteed to be written by the Rust
                   version of MLEM)
  --show-header    Show size headers included in self-describing raw 3D image
                   files, without visualizing the file in viewer
  --add-header     Add headers to headerless image files. The original files
                   will not be modified. The headers will be added to new
                   files matching the originals with the suffix
"""

import struct
from operator import itemgetter
from docopt import docopt

from utils import read_raw, read_raw_without_header, wrap_1d_into_3d, plt, np

class image:

    def __init__(self, data, shape, full_lengths):
        self.data = data
        self.shape = shape
        self.extents = full_lengths
        self.half_size = tuple(l/2 for l in full_lengths)
        self.pixel_size = tuple(l/n for l,n in zip(full_lengths, shape))


def image_and_voxel_size(voxel_values, n_voxels, fov_size):
    img  = image(wrap_1d_into_3d(voxel_values, n_voxels), n_voxels, fov_size)
    voxsize = tuple(s/n for n,s in zip(n_voxels, fov_size))
    return img, voxsize

class view:

    def __init__(self, files, *, axis='z', header=None, end='>'):
        self.files = files
        self.ts = set('z')
        self.fig, self.ax = plt.subplots()
        if header:
            shape, full_lengths = header
            def read_file(f): return read_raw_without_header(f, end=end)
            tmp = tuple(image_and_voxel_size(data, shape, full_lengths)
                        for                        data in map(read_file, files))
        else:
            def read_file(f): return read_raw(f, end=end)
            tmp = tuple(image_and_voxel_size(data, shape, full_lengths)
                        for (shape, full_lengths), data in map(read_file, files))

        self.images, self.voxel_sizes = tuple(zip(*tmp))
        self.maxima = [im.data.max() for im in self.images]
        self.image_number = 0
        self.axis = axis
        self.integrate = 0
        shape = self.images[self.image_number].shape
        self.pos = [n // 2 for n in shape] # start off in the middle along each axis
        self.ax.imshow(self.data_to_be_shown())
        self.aximage = self.ax.images[0]
        self.fig.canvas.mpl_connect('key_press_event', self.process_key)
        self.update()
        plt.show()

    def data_to_be_shown(self):
        s = [slice(None) for _ in range(3)]
        naxis = self.naxis
        n, i = self.pos[naxis], self.integrate
        s[naxis] = slice(n-i, n+i+1)
        image = self.images[self.image_number]
        it = image.data[s[0], s[1], s[2]]
        it = it.sum(axis = naxis)
        if self.axis in self.ts:
            it = it.transpose()
        it = np.flip(it,0)
        return it

    def process_key(self, event):
        k = event.key
        if k in ('p',       'down'): self.change_slice(-1)
        if k in ('P', 'shift+down'): self.change_slice(-10)
        if k in ('n',       'up')  : self.change_slice(+1)
        if k in ('N', 'shift+up')  : self.change_slice(+10)
        if k in 'xyz': self.axis = k
        if k == 't': self.flipT()
        if k in ('right', 'left'): self.switch_image(1 if k == 'right' else -1)
        if k == 'i': self.integrate = max(self.integrate - 1, 0)
        if k == 'I': self.integrate =     self.integrate + 1

        if k in 'c0': self.set_slice(0)
        self.update()

    def switch_image(self, n):
        self.image_number = (self.image_number + n) % len(self.images)
        self.update()

    def flipT(self):
        ax = self.axis
        if ax in self.ts: self.ts.remove(ax)
        else            : self.ts.add   (ax)

    def update(self):
        im = self.aximage
        data = self.data_to_be_shown()
        im.set_array(data)
        half_size = self.images[self.image_number].half_size
        ax, (x,y,z) = self.axis, half_size
        xe, ye,  xl, yl = ((y,x, 'y', 'x') if ax == 'z' else
                           (z,x, 'z', 'x') if ax == 'y' else
                           (z,y, 'z', 'y'))
        if self.axis in self.ts:
            xe,ye, xl,yl = ye,xe, yl,xl
        extent = -xe,xe, -ye,ye
        im.set_extent(extent)
        nax = self.naxis
        im.set_clim(0, data.max())
        #self.ax.set_aspect(ye/xe)
        i = self.pos[nax]
        pixel_size = self.images[self.image_number].pixel_size
        p = 'integrated' if self.integrate else f'{(i + 0.5) * pixel_size[nax] - half_size[nax]:6.1f}'
        nx,ny,nz = self.voxel_sizes[self.image_number]
        self.ax.set_title(f'''{self.axis} = {p}        voxel size = {nx:.2} x {ny:.2} x {nz:.2} mm
        {self.files[self.image_number]}''')
        self.ax.set_xlabel(xl)
        self.ax.set_ylabel(yl)
        self.fig.canvas.draw()

    def change_slice(self, n):
        ax = self.naxis
        self.pos[ax] = self.wrap(self.pos[ax] + n)

    def set_slice(self, n):
        ax = self.naxis
        self.pos[ax] = self.wrap(self.images[self.image_number].shape[ax] // 2 + n)

    def wrap(self, slice_):
        "Wrap around edges"
        return slice_ % self.images[self.image_number].shape[self.naxis]


    @property
    def naxis(self):
        ax = self.axis
        return 2 if ax == 'z' else 1 if ax == 'y' else 0


def add_header(filenames):
    for name in filenames:
        in_ = open(name     , 'rb').read()
        out = open(name+'_h', 'wb')
        pixels = struct.pack('>HHH', nx, ny, nz)
        mm     = struct.pack('>fff', dx, dy, dz)
        out.write(pixels)
        out.write(mm)
        out.write(in_)


if __name__ == '__main__':
    args = docopt(__doc__)

    endianness = '<' if args['--little-endian'] else '>'
    filenames = args['FILES']

    # Parse header from CLI, if necessary
    if args['--add-header'] or args['--assume-header']:
        nx,ny,nz = map(int  , itemgetter('NX','NY','NZ')(args))
        dx,dy,dz = map(float, itemgetter('DX','DY','DZ')(args))
        header_on_cli = (nx,ny,nz), (dx,dy,dz)
    else:
        header_on_cli = None

    if args['--add-header']:
        add_header(filenames)
        exit(0)

    if args['--show-header']:
        longest = len(max(filenames, key=len))
        for filename in filenames:
            (nx,ny,nz), (dx,dy,dz) = read_raw(filename, header_only=True)
            print(f'{filename:{longest}}:   {nx} {ny} {nz}   {dx} {dy} {dz}')
        exit(0)

    v = view(filenames, header=header_on_cli, end=endianness)

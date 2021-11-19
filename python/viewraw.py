import struct
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

    def __init__(self, files, *, axis='z', header=None):
        self.files = files
        self.ts = set('z')
        self.fig, self.ax = plt.subplots()
        if header:
            shape, full_lengths = header
            tmp = tuple(image_and_voxel_size(data, shape, full_lengths)
                        for                        data in map(read_raw_without_header, files))
        else:
            tmp = tuple(image_and_voxel_size(data, shape, full_lengths)
                        for (shape, full_lengths), data in map(read_raw,                files))

        self.images, self.voxel_sizes = tuple(zip(*tmp))
        self.maxima = [im.data.max() for im in self.images]
        self.image_number = 0
        self.axis = axis
        self.integrate = False
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
        if not self.integrate:
            s[naxis] = self.pos[naxis]
        image = self.images[self.image_number]
        it = image.data[s[0], s[1], s[2]]
        if self.integrate:
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
        if k == 'i': self.integrate = not self.integrate
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
        if self.integrate:
            im.set_clim(0, data.max())
        else:
            im.set_clim(0, self.   maxima[self.image_number])
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


def usage(argv):
    return f"""Usage:

    {argv[0]} FILENAMES...

       View raw 3D image files which contain a self-describing size header.

    {argv[0]} --show-header FILENAMES...

       Show size headers included in self-describing raw 3D image files.

    {argv[0]} --assume-header  NX NY NZ   DX DY DZ  FILENAMES...

       View headerless raw 3D image. Specify size on command line.

       NX, NY, NZ: number of voxels in each dimension
       dx, DY, DZ: full extent of FOV in each dimension

    {argv[0]} --add-header   NX NY NZ   DX DY DZ  FILENAMES...

       Add headers to headerless raw 3D image files.

       The original files will not be modified. The headers will be added to new
       files matching the original names with the suffix `_h`
    """


if __name__ == '__main__':
    from sys import argv

    if '-h' in argv or '--help' in argv:
        exit(usage(argv))
    filenames = list(f for f in argv[1:] if not f.startswith('-'))

    header_on_cli = None

    if argv[1] in ('--add-header', '--assume-header'):
        prog_, flag, nx,ny,nz, dx,dy,dz, *filenames = argv
        nx,ny,nz = map(int  , (nx,ny,nz))
        dx,dy,dz = map(float, (dx,dy,dz))

        if flag == '--add-header':
            add_header(filenames)
            exit(0)

        else:
            header_on_cli = (nx,ny,nz), (dx,dy,dz)

    if '--show-header' in argv:
        longest = len(max(filenames, key=len))
        for filename in filenames:
            (nx,ny,nz), (dx,dy,dz) = read_raw(filename, header_only=True)
            print(f'{filename:{longest}}:   {nx} {ny} {nz}   {dx} {dy} {dz}')

    else:
        v = view(filenames, header=header_on_cli)
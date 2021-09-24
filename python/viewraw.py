from utils import read_raw, wrap_1d_into_3d, plt, np

class image:

    def __init__(self, data, shape, full_lengths):
        self.data = data
        self.shape = shape
        self.extents = full_lengths
        self.half_size = tuple(l/2 for l in full_lengths)
        self.pixel_size = tuple(l/n for l,n in zip(full_lengths, shape))


class view:

    def __init__(self, files, *, axis='z'):
        self.files = files
        self.ts = set('z')
        self.fig, self.ax = plt.subplots()
        self.images = tuple(image(wrap_1d_into_3d(data, shape), shape, full_lengths)
                            for (shape, full_lengths), data in map(read_raw, files))

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
        self.ax.set_title(f'{self.axis} = {p}\n{self.files[self.image_number]}')
        self.ax.set_xlabel(xl)
        self.ax.set_ylabel(yl)
        self.fig.canvas.draw()

    def change_slice(self, n):
        ax = self.naxis
        self.pos[ax] = (self.pos[ax] + n) % self.images[self.image_number].shape[ax] # wrap around at edges

    @property
    def naxis(self):
        ax = self.axis
        return 2 if ax == 'z' else 1 if ax == 'y' else 0

if __name__ == '__main__':
    from sys import argv
    filenames = list(f for f in argv[1:] if not f.startswith('--'))
    if '--header-only' in argv:
        for filename in filenames:
            print(read_raw(filename, header_only=True))
    else:
        v = view(filenames)

from utils import read_raw, wrap_1d_into_3d, plt

class view:

    def __init__(self, *, data, shape, full_size, axis='x'):
        self.ts = set('z')
        self.fig, self.ax = plt.subplots()
        self.data3d = wrap_1d_into_3d(data, shape)
        self.shape = shape
        self.full_size = full_size
        self.half_size = tuple(l/2 for l in full_size)
        self.pixel_size = [full_size[n] / shape[n] for n in range(3)]
        self.axis = axis
        self.pos = [n // 2 for n in shape] # start off in the middle along each axis
        self.ax.imshow(self.slice_through_data())
        self.aximage = self.ax.images[0]
        self.fig.canvas.mpl_connect('key_press_event', self.process_key)
        self.update()
        plt.show()

    def slice_through_data(self):
        s = [slice(None) for _ in range(3)]
        naxis = self.naxis
        s[naxis] = self.pos[naxis]
        it = self.data3d[s[0], s[1], s[2]]
        if self.axis in self.ts:
            it = it.transpose()
        return it

    def process_key(self, event):
        k = event.key
        if k == 'p': self.change_slice(-1)
        if k == 'n': self.change_slice(+1)
        if k == 'P': self.change_slice(-10)
        if k == 'N': self.change_slice(+10)
        if k in 'xyz': self.axis = k
        if k == 't': self.flipT()
        self.update()

    def flipT(self):
        ax = self.axis
        if ax in self.ts: self.ts.remove(ax)
        else            : self.ts.add   (ax)

    def update(self):
        im = self.aximage
        im.set_array(self.slice_through_data())
        ax, (x,y,z) = self.axis, self.half_size
        xe, ye,  xl, yl = ((x,y, 'x', 'y') if ax == 'z' else
                           (z,x, 'z', 'x') if ax == 'y' else
                           (z,y, 'z', 'y'))
        extent = -xe,xe, -ye,ye
        im.set_extent(extent)
        #self.ax.set_aspect(ye/xe)
        nax = self.naxis
        i = self.pos[nax]
        p = (i + 0.5) * self.pixel_size[nax] - self.half_size[nax]
        self.ax.set_title(f'{self.axis} = {p}     T={"".join(sorted(self.ts))}')
        self.ax.set_xlabel(xl)
        self.ax.set_ylabel(yl)
        self.fig.canvas.draw()

    def change_slice(self, n):
        ax = self.naxis
        self.pos[ax] = (self.pos[ax] + n) % self.shape[ax] # wrap around at edges

    @property
    def naxis(self):
        ax = self.axis
        return 2 if ax == 'z' else 1 if ax == 'y' else 0

if __name__ == '__main__':
    from sys import argv

    (shape, full_length), data = read_raw(argv[1])
    view(data=data, shape=shape, full_size=full_length)

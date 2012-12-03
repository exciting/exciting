from ase.io.eps import EPS


class PNG(EPS):
    def write_header(self):
        from matplotlib.backends.backend_agg import RendererAgg

        try:
            from matplotlib.transforms import Value
        except ImportError:
            dpi = 72
        else:
            dpi = Value(72)

        self.renderer = RendererAgg(self.w, self.h, dpi)

        #self.gc = GraphicsContextBase()
        #self.gc.set_linewidth(2)

    def write_trailer(self):
        renderer = self.renderer
        if hasattr(renderer._renderer, 'write_png'):
            # Old version of matplotlib:
            renderer._renderer.write_png(self.filename)
        else:
            x = renderer._renderer.buffer_rgba(0, 0)
            from matplotlib import _png
            _png.write_png(renderer._renderer.buffer_rgba(0, 0),
                           renderer.width, renderer.height,
                           self.filename, 72)


def write_png(filename, atoms, **parameters):
    if isinstance(atoms, list):
        if len(atoms) > 1:
            raise RuntimeError("Don't know how to save more than "+
                               "one image to PNG image!")
        else:
            atoms = atoms[0]
    PNG(atoms, **parameters).write(filename)

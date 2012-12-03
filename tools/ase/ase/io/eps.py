import time
from math import sqrt

import numpy as np

from ase.utils import rotate
from ase.data import covalent_radii
from ase.data.colors import jmol_colors


class EPS:
    def __init__(self, atoms,
                 rotation='', show_unit_cell=False, radii=None,
                 bbox=None, colors=None, scale=20):
        self.numbers = atoms.get_atomic_numbers()
        self.colors = colors
        if colors is None:
            self.colors = jmol_colors[self.numbers]

        if radii is None:
            radii = covalent_radii[self.numbers]
        elif type(radii) is float:
            radii = covalent_radii[self.numbers] * radii
        else:
            radii = np.array(radii)

        natoms = len(atoms)

        if isinstance(rotation, str):
            rotation = rotate(rotation)

        A = atoms.get_cell()
        if show_unit_cell > 0:
            L, T, D = self.cell_to_lines(A)
            C = np.empty((2, 2, 2, 3))
            for c1 in range(2):
                for c2 in range(2):
                    for c3 in range(2):
                        C[c1, c2, c3] = np.dot([c1, c2, c3], A)
            C.shape = (8, 3)
            C = np.dot(C, rotation)  # Unit cell vertices
        else:
            L = np.empty((0, 3))
            T = None
            D = None
            C = None

        nlines = len(L)

        X = np.empty((natoms + nlines, 3))
        R = atoms.get_positions()
        X[:natoms] = R
        X[natoms:] = L

        r2 = radii**2
        for n in range(nlines):
            d = D[T[n]]
            if ((((R - L[n] - d)**2).sum(1) < r2) &
                (((R - L[n] + d)**2).sum(1) < r2)).any():
                T[n] = -1

        X = np.dot(X, rotation)
        R = X[:natoms]

        if bbox is None:
            X1 = (R - radii[:, None]).min(0)
            X2 = (R + radii[:, None]).max(0)
            if show_unit_cell == 2:
                X1 = np.minimum(X1, C.min(0))
                X2 = np.maximum(X2, C.max(0))
            M = (X1 + X2) / 2
            S = 1.05 * (X2 - X1)
            w = scale * S[0]
            if w > 500:
                w = 500
                scale = w / S[0]
            h = scale * S[1]
            offset = np.array([scale * M[0] - w / 2, scale * M[1] - h / 2, 0])
        else:
            w = (bbox[2] - bbox[0]) * scale
            h = (bbox[3] - bbox[1]) * scale
            offset = np.array([bbox[0], bbox[1], 0]) * scale

        self.w = w
        self.h = h

        X *= scale
        X -= offset

        if nlines > 0:
            D = np.dot(D, rotation)[:, :2] * scale

        if C is not None:
            C *= scale
            C -= offset

        A = np.dot(A, rotation)
        A *= scale

        self.A = A
        self.X = X
        self.D = D
        self.T = T
        self.C = C
        self.natoms = natoms
        self.d = 2 * scale * radii

    def cell_to_lines(self, A):
        nlines = 0
        nn = []
        for c in range(3):
            d = sqrt((A[c]**2).sum())
            n = max(2, int(d / 0.3))
            nn.append(n)
            nlines += 4 * n

        X = np.empty((nlines, 3))
        T = np.empty(nlines, int)
        D = np.zeros((3, 3))

        n1 = 0
        for c in range(3):
            n = nn[c]
            dd = A[c] / (4 * n - 2)
            D[c] = dd
            P = np.arange(1, 4 * n + 1, 4)[:, None] * dd
            T[n1:] = c
            for i, j in [(0, 0), (0, 1), (1, 0), (1, 1)]:
                n2 = n1 + n
                X[n1:n2] = P + i * A[(c + 1) % 3] + j * A[(c + 2) % 3]
                n1 = n2

        return X, T, D

    def write(self, filename):
        self.filename = filename
        self.write_header()
        self.write_body()
        self.write_trailer()

    def write_header(self):
        import matplotlib
        if matplotlib.__version__ <= '0.8':
            raise RuntimeError('Your version of matplotlib (%s) is too old' %
                               matplotlib.__version__)

        from matplotlib.backends.backend_ps import RendererPS, \
             GraphicsContextPS, psDefs

        self.fd = open(self.filename, 'w')
        self.fd.write('%!PS-Adobe-3.0 EPSF-3.0\n')
        self.fd.write('%%Creator: G2\n')
        self.fd.write('%%CreationDate: %s\n' % time.ctime(time.time()))
        self.fd.write('%%Orientation: portrait\n')
        bbox = (0, 0, self.w, self.h)
        self.fd.write('%%%%BoundingBox: %d %d %d %d\n' % bbox)
        self.fd.write('%%EndComments\n')

        Ndict = len(psDefs)
        self.fd.write('%%BeginProlog\n')
        self.fd.write('/mpldict %d dict def\n' % Ndict)
        self.fd.write('mpldict begin\n')
        for d in psDefs:
            d = d.strip()
            for l in d.split('\n'):
                self.fd.write(l.strip() + '\n')
        self.fd.write('%%EndProlog\n')

        self.fd.write('mpldict begin\n')
        self.fd.write('%d %d 0 0 clipbox\n' % (self.w, self.h))

        self.renderer = RendererPS(self.w, self.h, self.fd)

    def write_body(self):
        try:
            from matplotlib.path import Path
        except ImportError:
            Path = None
            from matplotlib.patches import Circle, Polygon
        else:
            from matplotlib.patches import Circle, PathPatch

        indices = self.X[:, 2].argsort()
        for a in indices:
            xy = self.X[a, :2]
            if a < self.natoms:
                r = self.d[a] / 2
                if ((xy[1] + r > 0) and (xy[1] - r < self.h) and
                    (xy[0] + r > 0) and (xy[0] - r < self.w)):
                    circle = Circle(xy, r, facecolor=self.colors[a])
                    circle.draw(self.renderer)
            else:
                a -= self.natoms
                c = self.T[a]
                if c != -1:
                    hxy = self.D[c]
                    if Path is None:
                        line = Polygon((xy + hxy, xy - hxy))
                    else:
                        line = PathPatch(Path((xy + hxy, xy - hxy)))
                    line.draw(self.renderer)

    def write_trailer(self):
        self.fd.write('end\n')
        self.fd.write('showpage\n')
        self.fd.close()


def write_eps(filename, atoms, **parameters):
    if isinstance(atoms, list):
        assert len(atoms) == 1
        atoms = atoms[0]
    EPS(atoms, **parameters).write(filename)

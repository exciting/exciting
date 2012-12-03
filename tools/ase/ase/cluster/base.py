import numpy as np

class ClusterBase:
    def get_layer_distance(self, miller, layers=1, tol=1e-9, new=True):
        """Returns the distance between planes defined by the given miller
        index.
        """
        if new:
            # Create lattice sample
            size = np.zeros(3, int)
            for i, m in enumerate(miller):
                size[i] = np.abs(m) + 2

            m = len(self.atomic_basis)
            p = np.zeros((size.prod() * m, 3))
            for h in range(size[0]):
                for k in range(size[1]):
                    for l in range(size[2]):
                        i = h * (size[1] * size[2]) + k * size[2] + l
                        p[m*i:m*(i+1)] = np.dot([h, k, l] + self.atomic_basis,
                                                self.lattice_basis)

            # Project lattice positions on the miller direction.
            n = self.miller_to_direction(miller)
            d = np.sum(n * p, axis=1)
            if np.all(d < tol):
                # All negative incl. zero
                d = np.sort(np.abs(d))
                reverse = True
            else:
                # Some or all positive
                d = np.sort(d[d > -tol])
                reverse = False
            d = d[np.concatenate((d[1:] - d[:-1] > tol, [True]))]
            d = d[1:] - d[:-1]

            # Look for a pattern in the distances between layers. A pattern is
            # accepted if more than 50 % of the distances obeys it.
            pattern = None
            for i in range(len(d)):
                for n in range(1, (len(d) - i) / 2 + 1):
                    if np.all(np.abs(d[i:i+n] - d[i+n:i+2*n]) < tol):
                        counts = 2
                        for j in range(i+2*n, len(d), n):
                            if np.all(np.abs(d[j:j+n] - d[i:i+n]) < tol):
                                counts += 1
                        if counts * n * 1.0 / len(d) > 0.5:
                            pattern = d[i:i+n].copy()
                            break
                if pattern is not None:
                    break

            if pattern is None:
                raise RuntimeError('Could not find layer distance for the ' +
                                   '(%i,%i,%i) surface.' % miller)
            if reverse:
                pattern = pattern[::-1]

            if layers < 0:
                pattern = -1 * pattern[::-1]
                layers *= -1

            map = np.arange(layers - layers % 1 + 1, dtype=int) % len(pattern)
            return pattern[map][:-1].sum() + layers % 1 * pattern[map][-1]



        n = self.miller_to_direction(miller)
        d1 = d2 = 0.0

        d = np.abs(np.sum(n * self.lattice_basis, axis=1))
        mask = np.greater(d, 1e-10)
        if mask.sum() > 0:
            d1 = np.min(d[mask])

        if len(self.atomic_basis) > 1:
            atomic_basis = np.dot(self.atomic_basis, self.lattice_basis)
            d = np.sum(n * atomic_basis, axis=1)
            s = np.sign(d)
            d = np.abs(d)
            mask = np.greater(d, 1e-10)
            if mask.sum() > 0:
                d2 = np.min(d[mask])
                s2 = s[mask][np.argmin(d[mask])]

        if d2 > 1e-10:
            if s2 < 0 and d1 - d2 > 1e-10:
                d2 = d1 - d2
            elif s2 < 0 and d2 - d1 > 1e-10:
                d2 = 2 * d1 - d2
            elif s2 > 0 and d2 - d1 > 1e-10:
                d2 = d2 - d1

            if np.abs(d1 - d2) < 1e-10:
                ld = np.array([d1])
            elif np.abs(d1 - 2 * d2) < 1e-10:
                ld = np.array([d2])
            else:
                assert d1 > d2, "Something is wrong with the layer distance."
                ld = np.array([d2, d1 - d2])
        else:
            ld = np.array([d1])

        if len(ld) > 1:
            if layers < 0:
                ld = np.array([-ld[1], -ld[0]])
                layers *= -1

            map = np.arange(layers - (layers % 1), dtype=int) % len(ld)
            r = ld[map].sum() + (layers % 1) * ld[np.abs(map[-1] - 1)]
        else:
            r = ld[0] * layers

        return r

    def miller_to_direction(self, miller, norm=True):
        """Returns the direction corresponding to a given Miller index."""
        d = np.dot(miller, self.resiproc_basis)
        if norm:
            d = d / np.linalg.norm(d)
        return d


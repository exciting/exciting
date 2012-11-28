from math import exp, pi, sin, sqrt, cos, acos
import numpy as np

from ase.data import atomic_numbers

# Table (1) of
# D. WAASMAIER AND A. KIRFEL, Acta Cryst. (1995). A51, 416-431
waasmaier = {
    #      a1        b1         a2        b2        a3        b3          a4         b4         a5         b5        c
    'C' : [2.657506, 14.780758, 1.078079, 0.776775, 1.490909, 42.086843,  -4.241070, -0.000294, 0.713791, 0.239535, 4.297983],
    'S' : [6.372157, 1.514347, 5.154568, 22.092528, 1.473732, 0.061373,   1.635073,  55.445176, 1.209372, 0.646925, 0.154722],
    'Pd': [6.121511, 0.062549,  4.784063, 0.784031, 16.631683, 8.751391,  4.318258, 34.489983, 13.246773, 0.784031, 0.883099],
    'Ag': [6.073874, 0.055333, 17.155437, 7.896512, 4.173344, 28.443739,  0.852238, 110.376108, 17.988685, 0.716809, 0.756603],
    'Au': [16.777389, 0.122737, 19.317156, 8.621570, 32.979682, 1.256902, 5.595453, 38.008821, 10.576854, 0.000601, -6.279078],
    'P' : [1.950541, 0.908139, 4.146930, 27.044953, 1.494560, 0.071280, 1.522042, 67.520190, 5.729711, 1.981173, 0.155233],
    'Cl': [1.446071, 0.052357, 6.870609, 1.193165, 6.151801, 18.343416, 1.750347, 46.398394, 0.634168, 0.401005, 0.146773],
}

class XrDebye:
    def __init__(self, wavelength, alpha=1.01, damping=0.04, warn=True,
                 method='Iwasa'):
        """
        Obtain powder x-ray spectra.

        wavelength in Angstrom
        damping in Angstrom**2
        """
        self.wavelength = wavelength
        self.damping = damping
        self.alpha = alpha
        self.warn = warn
        self.method = method

    def set_damping(self, damping):
        self.damping = damping

    def get(self, atoms, s):
        """Get the powder x-ray (XRD) pattern using the Debye-Formula.

        After: T. Iwasa and K. Nobusada, J. Phys. Chem. C 111 (2007) 45
               s is assumed to be in 1/Angstrom
        """

        sinth = self.wavelength * s / 2.
        costh = sqrt(1. - sinth**2)
        cos2th = cos(2. * acos(costh))
        pre = exp(- self.damping * s**2 / 2)
 
        if self.method == 'Iwasa':
            pre *= costh / (1. + self.alpha * cos2th**2)

        f = {}
        def atomic(symbol):
            if not f.has_key(symbol):
                if self.method == 'Iwasa':
                    f[symbol] = self.get_waasmaier(symbol, s)
                else:
                    f[symbol] = atomic_numbers[symbol]
            return f[symbol]

        def sinc(x):
            if x < 1.e-6:
                x2 = x * x
                return 1 - x2 / 6. + x2 * x2 / 120.
            else:
                return sin(x) / x

        I = 0.
        for a in atoms:
            fa = atomic(a.symbol)
#            print a.symbol, fa
            for b in atoms:
                fb = atomic(b.symbol)

                if a == b:
                    twopis = 0.
                else:
                    vrij = a.position - b.position
                    rij = np.sqrt(np.dot(vrij, vrij))
                    twopisr = 2 * pi * s * rij

                I += fa * fb * sinc(twopisr)
                    
        return pre * I

    def get_waasmaier(self, symbol, s):
        """Scattering factor for free atoms."""
        if symbol == 'H':
            # XXXX implement analytical H
            return 0
        elif waasmaier.has_key(symbol):
            abc = waasmaier[symbol]
            f = abc[10]
            s2 = s*s
            for i in range(5):
                f += abc[2 * i] * exp(-abc[2 * i + 1] * s2)
            return f
        if self.warn:
            print '<xrdebye::get_atomic> Element', symbol, 'not available'
        return 0
        

from math import pi, sqrt

# Constants from Konrad Hinsen's PhysicalQuantities module (1986 CODATA):
_c = 299792458.              # speed of light, m/s
_mu0 = 4.e-7 * pi            # permeability of vacuum
_eps0 = 1 / _mu0 / _c**2     # permittivity of vacuum
_Grav = 6.67259e-11          # gravitational constant
_hplanck = 6.6260755e-34     # Planck constant, J s
_hbar = _hplanck / (2 * pi)  # Planck constant / 2pi, J s
_e = 1.60217733e-19          # elementary charge
_me = 9.1093897e-31          # electron mass
_mp = 1.6726231e-27          # proton mass
_Nav = 6.0221367e23          # Avogadro number
_k = 1.380658e-23            # Boltzmann constant, J/K
_amu = 1.6605402e-27         # atomic mass unit, kg

Ang = Angstrom = 1.0
nm = 10.0
Bohr = 4e10 * pi * _eps0 * _hbar**2 / _me / _e**2  # Bohr radius

eV = 1.0
Hartree = _me * _e**3 / 16 / pi**2 / _eps0**2 / _hbar**2
kJ = 1000.0 / _e
kcal = 4.184 * kJ
mol = _Nav
Rydberg = 0.5 * Hartree
Ry = Rydberg
Ha = Hartree

second = 1e10 * sqrt(_e / _amu)
fs = 1e-15 * second

kB = _k / _e                 # Boltzmann constant, eV/K

Pascal = (1 / _e) / 1e30  # J/m^3
GPa = 1e9 * Pascal

Debye = 1.0 / 1e11 / _e / _c
alpha = _e**2 / (4 * pi * _eps0) / _hbar / _c # fine structure constant

# Derived atomic units that have no assigned name:
_aut = _hbar / (alpha**2 * _me * _c**2)      # atomic unit of time, s
_auv =  _e**2 / _hbar / (4 * pi * _eps0)     # atomic unit of velocity, m/s
_auf = alpha**3 * _me**2 * _c**3 / _hbar     # atomic unit of force, N
_aup = alpha**5 * _me**4 * _c**5 / _hbar**3  # atomic unit of pressure, Pa

AUT = second * _aut

# SI units
m = 1e10 * Ang    # metre
kg = 1. / _amu    # kilogram
s = second        # second
A = 1.0 / _e / s  # ampere
# derived
J = kJ / 1000  # Joule   = kg * m**2 / s**2
C = 1.0 / _e   # Coulomb = A * s

del pi, sqrt

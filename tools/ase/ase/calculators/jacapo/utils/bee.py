'''
adapted from ASE/Utilities/BEE.py
'''

from math import sqrt
import numpy as num
import numpy.random as ra

"""Bayesian Error Estimation

For details, see: "Bayesian Error Estimation in Density Functional
Theory", J. J. Mortensen, K. Kaasbjerg, S. L. Frederiksen,
J. K. Norskov, J. P. Sethna, K. W. Jacobsen, Phys. Rev. Lett. 95,
216401 (2005)."""

#                                 T
# cost(c) = cost0 + 0.5 * (c - c0) H (c - c0)
#

# Cost function minimum value:
cost0 = 3.4660625596

# Best fit parameters:
c0 = num.array([1.000787451, 0.1926284063, 1.896191546])

# Hessian:
# H = num.array([[ 1.770035168e+03, -3.732470432e+02, -2.105836167e+02],
#                [-3.732470432e+02,  1.188857209e+02,  6.054102443e+01],
#                [-2.105836167e+02,  6.054102443e+01,  3.211200293e+01]])
#
# 0.5 * np * T = cost0 (np=3: number of parameters)
T = cost0 * 2 / 3

def MakeEnsemble(N=1000, seed=117):
    ra.seed(seed)
    M = num.array([(0.066, -0.812, 1.996),
                   (0.055, 0.206, 0.082),
                   (-0.034, 0.007, 0.004)])
    alpha = ra.normal(0.0, 1.0, (N, 3))
    return c0 + num.dot(alpha, M)

c = MakeEnsemble()

def GetEnsembleEnergies(atoms, c=c):
    if hasattr(atoms, 'get_calculator'):
        coefs = atoms.get_calculator().get_ensemble_coefficients()
    else:
        coefs = atoms
    return coefs[0] + num.dot(c, coefs[1:])

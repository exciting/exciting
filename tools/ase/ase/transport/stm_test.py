import numpy as np
import pylab
import ase.transport.stm as stm


#           Parameters for a simple model.
#
#                           
#                         * eps_a
#                 v_ts /   \ v_a2 
#  ... *  *  *  *            *  *  *  *  * ...
#       \/                          \/
#       t1                          t2
#
#       Tip                      Surface
# ----------------|      |-----------------------
t1 = -1.0
t2 = -2.0
eps_a = 0.4

v_ts  = 0.05
v_a2 = 1.0

#Tip
h1 = np.zeros([2, 2])
h1[0, 1] = t1
h1[1, 0] = t1
s1 = np.identity(2)

h10 = np.zeros([2,2])
h10[0, 1] = t1
h10[1, 0] = t1
s10 = np.identity(2)


#Surface with "molecule" a.
h2 = np.zeros([2,2])
h2[0, 1] = v_a2 
h2[1, 0] = v_a2
h1[0, 0] = eps_a
s2 = np.identity(2)

h20 = np.zeros([2,2])
h20[0, 1] = t2
h20[1, 0] = t2
s20 = np.identity(2)


#Tip Surface coupling
V_ts = np.zeros([2,2])
V_ts[1, 0] = v_ts

eta1 = 0.0001
eta2 = 0.0001

stm = reload(stm)
stm_calc = stm.STM(h1, s1, h2, s2, h10, s10, h20, s20, eta1, eta2)
energies = np.arange(-3.0, 3.0, 0.01)
stm_calc.initialize(energies)


T_stm = stm_calc.get_transmission(V_ts)


#Perform the full calculation and compare
from ase.transport.calculators import TransportCalculator as TC

h = np.zeros([4,4])
h[:2, :2] = h1
h[-2:, -2:] = h2
h[:2, -2:] = V_ts
h[-2:, :2] = V_ts.T


tc = TC(energies=energies,
        h=h,
        h1=h10,
        h2=h20,
        eta=eta1, eta1=eta1, eta2=eta2)
        
T_full = tc.get_transmission()

pylab.plot(stm_calc.energies, T_stm, 'b')
pylab.plot(tc.energies, T_full, 'r--')
pylab.show()


#bias stuff
biass = np.arange(-2.0, 2.0, 0.2)
Is = [stm_calc.get_current(bias, V_ts) for bias in biass]
pylab.plot(biass, Is, '+')
pylab.show()



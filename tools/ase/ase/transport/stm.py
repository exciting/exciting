import time

import numpy as np

from ase.transport.tools import dagger
from ase.transport.selfenergy import LeadSelfEnergy
from ase.transport.greenfunction import GreenFunction
from ase.parallel import world


class STM:
    def __init__(self, h1, s1, h2, s2 ,h10, s10, h20, s20, eta1, eta2, w=0.5, pdos=[], logfile = None):
        """XXX
        
        1. Tip
        2. Surface
        
        h1: ndarray
            Hamiltonian and overlap matrix for the isolated tip
            calculation.  Note, h1 should contain (at least) one
            principal layer.

        h2: ndarray
            Same as h1 but for the surface.

        h10: ndarray
            periodic part of the tip. must include two and only
            two principal layers.

        h20: ndarray
            same as h10, but for the surface

        The s* are the corresponding overlap matrices.  eta1, and eta
        2 are (finite) infinitesimals.  """
        
        self.pl1 = len(h10) // 2 #principal layer size for the tip
        self.pl2 = len(h20) // 2 #principal layer size for the surface
        self.h1 = h1 
        self.s1 = s1
        self.h2 = h2
        self.s2 = s2
        self.h10 = h10 
        self.s10 = s10 
        self.h20 = h20 
        self.s20 = s20
        self.eta1 = eta1
        self.eta2 = eta2
        self.w = w #asymmetry of the applied bias (0.5=>symmetric)
        self.pdos = []
        self.log = logfile

    def initialize(self, energies, bias=0):
        """
            energies: list of energies 
            for which the transmission function should be evaluated.
            bias.
            Will precalculate the surface greenfunctions of the tip and
            surface.
        """
        self.bias = bias
        self.energies = energies
        nenergies = len(energies)
        pl1, pl2 = self.pl1, self.pl2
        nbf1, nbf2 = len(self.h1), len(self.h2)
      
        #periodic part of the tip
        hs1_dii = self.h10[:pl1, :pl1], self.s10[:pl1, :pl1]
        hs1_dij = self.h10[:pl1, pl1:2*pl1], self.s10[:pl1, pl1:2*pl1]
        #coupling betwen per. and non. per part of the tip
        h1_im = np.zeros((pl1, nbf1), complex) 
        s1_im = np.zeros((pl1, nbf1), complex)
        h1_im[:pl1, :pl1], s1_im[:pl1, :pl1] = hs1_dij
        hs1_dim = [h1_im, s1_im]

        #periodic part the surface 
        hs2_dii = self.h20[:pl2, :pl2], self.s20[:pl2, :pl2]
        hs2_dij = self.h20[pl2:2*pl2, :pl2], self.s20[pl2:2*pl2, :pl2]
        #coupling betwen per. and non. per part of the surface
        h2_im = np.zeros((pl2, nbf2), complex)
        s2_im = np.zeros((pl2, nbf2), complex) 
        h2_im[-pl2:, -pl2:], s2_im[-pl2:, -pl2:] = hs2_dij
        hs2_dim = [h2_im, s2_im]

        #tip and surface greenfunction 
        self.selfenergy1 = LeadSelfEnergy(hs1_dii, hs1_dij, hs1_dim, self.eta1)
        self.selfenergy2 = LeadSelfEnergy(hs2_dii, hs2_dij, hs2_dim, self.eta2)
        self.greenfunction1 = GreenFunction(self.h1-self.bias*self.w*self.s1, self.s1, 
                                            [self.selfenergy1], self.eta1)
        self.greenfunction2 = GreenFunction(self.h2-self.bias*(self.w-1)*self.s2, self.s2, 
                                            [self.selfenergy2], self.eta2)
        
        #Shift the bands due to the bias.
        bias_shift1 = -bias * self.w
        bias_shift2 = -bias * (self.w - 1)
        self.selfenergy1.set_bias(bias_shift1)
        self.selfenergy2.set_bias(bias_shift2)
        
        #tip and surface greenfunction matrices.
        nbf1_small = nbf1 #XXX Change this for efficiency in the future
        nbf2_small = nbf2 #XXX -||-
        coupling_list1 = range(nbf1_small)# XXX -||-
        coupling_list2 = range(nbf2_small)# XXX -||-
        self.gft1_emm = np.zeros((nenergies, nbf1_small, nbf1_small), complex) 
        self.gft2_emm = np.zeros((nenergies, nbf2_small, nbf2_small), complex)
 
        for e, energy in enumerate(self.energies):
            if self.log != None: # and world.rank == 0:
                    T = time.localtime()
                    self.log.write(' %d:%02d:%02d, ' % (T[3], T[4], T[5]) +
                                   '%d, %d, %02f\n' % (world.rank, e, energy))
            gft1_mm = self.greenfunction1.retarded(energy)[coupling_list1]
            gft1_mm = np.take(gft1_mm, coupling_list1, axis=1)

            gft2_mm = self.greenfunction2.retarded(energy)[coupling_list2]
            gft2_mm = np.take(gft2_mm, coupling_list2, axis=1)
 
            self.gft1_emm[e] = gft1_mm
            self.gft2_emm[e] = gft2_mm

            if self.log != None and world.rank == 0:
                self.log.flush()

    def get_transmission(self, v_12, v_11_2=None, v_22_1=None):
        """XXX

        v_12:
            coupling between tip and surface 
        v_11_2:
            correction to "on-site" tip elements due to the 
            surface (eq.16). Is only included to first order.
        v_22_1:
            corretion to "on-site" surface elements due to he
            tip (eq.17). Is only included to first order.
        """

        dim0 = v_12.shape[0]
        dim1 = v_12.shape[1]

        nenergies = len(self.energies)
        T_e = np.empty(nenergies,float)
        v_21 = dagger(v_12)
        for e, energy in enumerate(self.energies):
            gft1 = self.gft1_emm[e]
            if v_11_2!=None:
                gf1 = np.dot(v_11_2, np.dot(gft1, v_11_2)) 
                gf1 += gft1 #eq. 16
            else:
                gf1 = gft1
            
            gft2 = self.gft2_emm[e]
            if v_22_1!=None:
                gf2 = np.dot(v_22_1,np.dot(gft2, v_22_1))
                gf2 += gft2 #eq. 17
            else:
                gf2 = gft2
            
            a1 = (gf1 - dagger(gf1))
            a2 = (gf2 - dagger(gf2))
            self.v_12 = v_12
            self.a2 = a2
            self.v_21 = v_21
            self.a1 = a1
            v12_a2 = np.dot(v_12, a2[:dim1])
            v21_a1 = np.dot(v_21, a1[-dim0:])
            self.v12_a2 = v12_a2
            self.v21_a1 = v21_a1
            T = -np.trace(np.dot(v12_a2[:,:dim1], v21_a1[:,-dim0:])) #eq. 11
            assert abs(T.imag).max() < 1e-14
            T_e[e] = T.real
            self.T_e = T_e
        return T_e


    def get_current(self, bias, v_12, v_11_2=None, v_22_1=None):
        """Very simple function to calculate the current.
        
        Asummes zero temperature.

        bias: type? XXX
            bias voltage (V)
            
        v_12: XXX
            coupling between tip and surface.
            
        v_11_2:
            correction to onsite elements of the tip
            due to the potential of the surface.
        v_22_1:
            correction to onsite elements of the surface
            due to the potential of the tip.
        """
        energies = self.energies
        T_e = self.get_transmission(v_12, v_11_2, v_22_1)
        bias_window = -np.array([bias * self.w, bias * (self.w - 1)])
        bias_window.sort()
        self.bias_window = bias_window
        #print 'bias window', np.around(bias_window,3)
        #print 'Shift of tip lead do to the bias:', self.selfenergy1.bias
        #print 'Shift of surface lead do to the bias:', self.selfenergy2.bias
        i1 = sum(energies < bias_window[0]) 
        i2 = sum(energies < bias_window[1])
        step = 1 
        if i2 < i1:
            step = -1
        
        return np.sign(bias)*np.trapz(x=energies[i1:i2:step], y=T_e[i1:i2:step])







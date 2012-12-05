# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.

"""
Quasi-Newton algorithm
"""

__docformat__ = 'reStructuredText'

import numpy as np
import weakref,time,sys


def f(lamda,Gbar,b,radius): 
        b1 = b - lamda
        g = radius**2 - np.dot(Gbar/b1, Gbar/b1)
        return g



def scale_radius_energy(f,r):
        scale = 1.0
#       if(r<=0.01):
#               return scale
        
        if f<0.01: scale*=1.4
        if f<0.05: scale*=1.4
        if f<0.10: scale*=1.4
        if f<0.40: scale*=1.4

        if f>0.5: scale *= 1./1.4               
        if f>0.7: scale *= 1./1.4               
        if f>1.0: scale *= 1./1.4

        return scale

def scale_radius_force(f,r):
        scale = 1.0
#       if(r<=0.01):
#               return scale
        g = abs(f -1)
        if g<0.01: scale*=1.4
        if g<0.05: scale*=1.4
        if g<0.10: scale*=1.4
        if g<0.40: scale*=1.4

        if g>0.5: scale *= 1./1.4               
        if g>0.7: scale *= 1./1.4               
        if g>1.0: scale *= 1./1.4

        return scale

def find_lamda(upperlimit,Gbar,b,radius):
        lowerlimit = upperlimit
        eps = 1e-12
        step = 0.1
        while  f(lowerlimit,Gbar,b,radius) < 0:
                lowerlimit -= step
                
        converged = False

        while not converged: 

                midt = (upperlimit+lowerlimit)/2.
                lamda = midt
                fmidt = f(midt,Gbar,b,radius)
                fupper = f(upperlimit,Gbar,b,radius)
                flower = f(lowerlimit,Gbar,b,radius)
        
                if fupper*fmidt<0: 
                        lowerlimit = midt 
                else: 
                        upperlimit = midt

                if abs(upperlimit-lowerlimit)<1e-6: 
                        converged = True

        return lamda

def get_hessian_inertia(eigenvalues):
        # return number of negative modes
        n = 0
        print 'eigenvalues ',eigenvalues[0],eigenvalues[1],eigenvalues[2]
        while eigenvalues[n]<0:
                n+=1
        return n 


from numpy.linalg import eigh, solve

from ase.optimize.optimize import Optimizer



class GoodOldQuasiNewton(Optimizer):

    def __init__(self, atoms, restart=None, logfile='-', trajectory=None,
                 fmax=None, converged=None,
                hessianupdate='BFGS',hessian=None,forcemin=True,
                verbosity=None,maxradius=None,
                diagonal=20.,radius=None,
                transitionstate = False):
            
        Optimizer.__init__(self, atoms, restart, logfile, trajectory)

        self.eps = 1e-12
        self.hessianupdate = hessianupdate
        self.forcemin = forcemin
        self.verbosity = verbosity
        self.diagonal = diagonal

        self.atoms = atoms

        n = len(self.atoms) * 3
        if radius is None: 
                self.radius = 0.05*np.sqrt(n)/10.0
        else:
                self.radius = radius

        if maxradius is None: 
                self.maxradius = 0.5*np.sqrt(n)
        else:
                self.maxradius = maxradius
                
        # 0.01 < radius < maxradius
        self.radius = max(min( self.radius, self.maxradius ), 0.0001)

        self.transitionstate = transitionstate

        # check if this is a nudged elastic band calculation
        if hasattr(atoms,'springconstant'): 
                self.forcemin=False

        self.t0 = time.time() 

    def initialize(self):pass

    def write_log(self,text):
        if self.logfile is not None:
            self.logfile.write(text + '\n')
            self.logfile.flush()

    def set_max_radius(self, maxradius):
                self.maxradius = maxradius
                self.radius = min(self.maxradius, self.radius)
                
    def set_hessian(self,hessian):
        self.hessian = hessian

    def get_hessian(self):
        if not hasattr(self,'hessian'): 
                self.set_default_hessian() 
        return self.hessian

    def set_default_hessian(self): 
        # set unit matrix
        n = len(self.atoms) * 3
        hessian = np.zeros((n,n)) 
        for i in range(n): 
                        hessian[i][i] = self.diagonal
        self.set_hessian(hessian) 

    def read_hessian(self,filename): 
        import cPickle
        f = open(filename,'r')
        self.set_hessian(cPickle.load(f))
        f.close()

    def write_hessian(self,filename): 
        import cPickle
        f = paropen(filename,'w')
        cPickle.dump(self.get_hessian(),f)
        f.close()

    def write_to_restartfile(self):
        import cPickle
        f = paropen(self.restartfile,'w')
        cPickle.dump((self.oldpos,
                      self.oldG,
                      self.oldenergy,
                      self.radius,
                      self.hessian,
                      self.energy_estimate),f)
        f.close()
        


    def update_hessian(self,pos,G):
        import copy
        if hasattr(self,'oldG'): 
                if self.hessianupdate=='BFGS': 
                        self.update_hessian_bfgs(pos,G) 
                elif self.hessianupdate== 'Powell': 
                        self.update_hessian_powell(pos,G) 
                else:           
                        self.update_hessian_bofill(pos,G) 
        else: 
                if not hasattr(self,'hessian'): 
                        self.set_default_hessian()

        self.oldpos = copy.copy(pos)
        self.oldG = copy.copy(G)

        if self.verbosity: 
                print 'hessian ',self.hessian


        
    def update_hessian_bfgs(self,pos,G): 
        n = len(self.hessian)
        dgrad = G - self.oldG
        dpos  = pos - self.oldpos
        absdpos = np.sqrt(np.dot(dpos, dpos))
        dotg  = np.dot(dgrad,dpos) 
        tvec  = np.dot(dpos,self.hessian)
        dott  = np.dot(dpos,tvec)
        if (abs(dott)>self.eps) and (abs(dotg)>self.eps): 
                for i in range(n): 
                        for j in range(n): 
                                h = dgrad[i]*dgrad[j]/dotg - tvec[i]*tvec[j]/dott
                                self.hessian[i][j] += h



    def update_hessian_powell(self,pos,G):          
        n = len(self.hessian)
        dgrad = G - self.oldG
        dpos  = pos - self.oldpos
        absdpos = np.dot(dpos, dpos)
        if absdpos<self.eps: 
                return

        dotg  = np.dot(dgrad,dpos) 
        tvec  = dgrad-np.dot(dpos,self.hessian)
        tvecdot = np.dot(tvec,tvec)
        tvecdpos = np.dot(tvec,dpos) 
        ddot = tvecdpos/absdpos

        dott  = np.dot(dpos,tvec)
        if (abs(dott)>self.eps) and (abs(dotg)>self.eps): 
                for i in range(n): 
                        for j in range(n): 
                                h = tvec[i]*dpos[j] + dpos[i]*tvec[j]-ddot*dpos[i]*dpos[j]
                                h *= 1./absdpos
                                self.hessian[i][j] += h


    def update_hessian_bofill(self,pos,G):                                                                     
        print 'update Bofill'
        n = len(self.hessian)                                                                               
        dgrad = G - self.oldG                                                                               
        dpos  = pos - self.oldpos                                                                           
        absdpos = np.dot(dpos, dpos)                                                                          
        if absdpos<self.eps: 
                return
        dotg  = np.dot(dgrad,dpos)                                                                         
        tvec  = dgrad-np.dot(dpos,self.hessian)                                                 
        tvecdot = np.dot(tvec,tvec)                                                                        
        tvecdpos = np.dot(tvec,dpos)                                                                       
        ddot = tvecdpos/absdpos                                                                             

        coef1 = 1. - tvecdpos*tvecdpos/(absdpos*tvecdot)
        coef2 = (1. - coef1)*absdpos/tvecdpos
        coef3 = coef1*tvecdpos/absdpos

        dott  = np.dot(dpos,tvec)                                                                          
        if (abs(dott)>self.eps) and (abs(dotg)>self.eps):                                                   
                for i in range(n):                                                                          
                        for j in range(n):                                                                  
                                h = coef1*(tvec[i]*dpos[j] + dpos[i]*tvec[j])-dpos[i]*dpos[j]*coef3 + coef2*tvec[i]*tvec[j]
                                h *= 1./absdpos
                                self.hessian[i][j] += h                                                     



    def step(self, f):
        """ Do one QN step
        """

        pos = self.atoms.get_positions().ravel()
        G = -self.atoms.get_forces().ravel()
        energy = self.atoms.get_potential_energy()


        self.write_iteration(energy,G)

        if hasattr(self,'oldenergy'):

                self.write_log('energies ' + str(energy) + ' ' + str(self.oldenergy))

                if self.forcemin:
                        de = 1e-4
                else:
                        de = 1e-2

                if self.transitionstate:
                        de = 0.2

                if (energy-self.oldenergy)>de:
                        self.write_log('reject step')
                        self.atoms.set_positions(self.oldpos.reshape((-1, 3)))
                        G = self.oldG
                        energy = self.oldenergy
                        self.radius *= 0.5
                else: 
                        self.update_hessian(pos,G)
                        de = energy - self.oldenergy
                        f = 1.0
                        if self.forcemin: 
                                self.write_log("energy change; actual: %f estimated: %f "%(de,self.energy_estimate))
                                if abs(self.energy_estimate)>self.eps: 
                                        f = abs((de/self.energy_estimate)-1)
                                        self.write_log('Energy prediction factor ' + str(f))
                                        # fg = self.get_force_prediction(G)
                                        self.radius *= scale_radius_energy(f,self.radius) 

                        else:
                                self.write_log("energy change; actual: %f "%(de))
                                self.radius*=1.5

                        fg = self.get_force_prediction(G)
                        self.write_log("Scale factors %f %f "%(scale_radius_energy(f,self.radius),
                                                                scale_radius_force(fg,self.radius)))
                        
                                   
                self.radius = max(min(self.radius,self.maxradius), 0.0001)
        else: 
                self.update_hessian(pos,G)

        self.write_log("new radius %f "%(self.radius))          
        self.oldenergy = energy

        b,V = eigh(self.hessian)
        V=V.T.copy()
        self.V = V

        # calculate projection of G onto eigenvectors V
        Gbar = np.dot(G,np.transpose(V))
        
        lamdas = self.get_lambdas(b,Gbar)

        D = -Gbar/(b-lamdas) 
        n = len(D)
        step = np.zeros((n))
        for i in range(n): 
                step += D[i]*V[i]

        pos = self.atoms.get_positions().ravel()
        pos += step

        energy_estimate = self.get_energy_estimate(D,Gbar,b) 
        self.energy_estimate = energy_estimate
        self.gbar_estimate = self.get_gbar_estimate(D,Gbar,b)
        self.old_gbar = Gbar

        self.atoms.set_positions(pos.reshape((-1, 3)))




    def get_energy_estimate(self,D,Gbar,b): 

        de = 0.0
        for n in range(len(D)): 
                de += D[n]*Gbar[n] + 0.5*D[n]*b[n]*D[n]
        return de

    def get_gbar_estimate(self,D,Gbar,b):
        gbar_est = (D*b) + Gbar
        self.write_log('Abs Gbar estimate ' + str(np.dot(gbar_est,gbar_est)))
        return gbar_est

    def get_lambdas(self,b,Gbar):
        lamdas = np.zeros((len(b)))

        D = -Gbar/b
        #absD = np.sqrt(np.sum(D**2))
        absD = np.sqrt(np.dot(D, D))

        eps = 1e-12
        nminus = self.get_hessian_inertia(b)

        if absD < self.radius:
                if not self.transitionstate:
                        self.write_log('Newton step') 
                        return lamdas
                else:
                        if nminus==1:
                                self.write_log('Newton step')
                                return lamdas
                        else:
                                self.write_log("Wrong inertia of Hessian matrix: %2.2f %2.2f "%(b[0],b[1]))

        else:
                self.write_log("Corrected Newton step: abs(D) = %2.2f "%(absD))

        if not self.transitionstate: 
                # upper limit
                upperlimit = min(0,b[0])-eps
                lowerlimit = upperlimit
                lamda = find_lamda(upperlimit,Gbar,b,self.radius)
                lamdas += lamda
        else:
                # upperlimit
                upperlimit = min(-b[0],b[1],0)-eps
                lamda = find_lamda(upperlimit,Gbar,b,self.radius)
                lamdas += lamda
                lamdas[0] -= 2*lamda
                
        return lamdas



    def print_hessian(self): 
        hessian = self.get_hessian()
        n = len(hessian)
        for i in range(n): 
            for j in range(n): 
                print "%2.4f " %(hessian[i][j]),
            print " "


    

    def get_hessian_inertia(self,eigenvalues):
        # return number of negative modes
        self.write_log("eigenvalues %2.2f %2.2f %2.2f "%(eigenvalues[0],
                                                        eigenvalues[1],
                                                        eigenvalues[2]))
        n = 0
        while eigenvalues[n]<0:
                n+=1
        return n

    def get_force_prediction(self,G):
        # return measure of how well the forces are predicted
        Gbar = np.dot(G,np.transpose(self.V))
        dGbar_actual = Gbar-self.old_gbar
        dGbar_predicted = Gbar-self.gbar_estimate

        f = np.dot(dGbar_actual,dGbar_predicted)/np.dot(dGbar_actual,dGbar_actual)
        self.write_log('Force prediction factor ' + str(f))
        return f

    def write_iteration(self,energy,G):pass

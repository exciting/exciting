import os
import numpy as np
import matplotlib.pyplot as plt
from ase.calculators.jacapo import *
from ase.dft.dos import DOS

class BandStructure:
    '''outline of class to facilitate band structure calculations
    '''
    def __init__(self,
                 atoms,
                 BZpath=[],
                 npoints=10,
                 outnc='harris.nc'):
        """Headline here ... XXX.
        
        atoms is an ase.Atoms object with calculator
        attached. Presumably the self-consistent charge density has
        already been calculated, otherwise, it will be.

        BZpath is a list of tuples describing the path through the
        Brillouin zone. The tuples have the form (label, kpt), e.g. ::

          [('$\Gamma$',[0.0, 0.0, 0.0]),
           ('X',[0.0, 0.5, 0.5]),
           ('L',[0.5, 0.0, 0.0]),
           ('$\Gamma$',[0.0, 0.0, 0.0])]

        the label is used in the figure and can include latex markup.

        npoints is the number of points on each segment. It can either
        be a constant, which is used for every segment, or a list of
        integers that is an integer for each segment.        
        """

        self.atoms = atoms
        self.calc = atoms.get_calculator()
        #first, we make sure the charge density is up to date.
        self.calc.get_charge_density()
        self.ef = self.calc.get_ef() #self-consistent fermi level

        self.labels = [x[0] for x in BZpath]
        self.kpt_path = [np.array(x[1],dtype=np.float) for x in BZpath]
        self.npoints = npoints

        #first, setup the kpt path
        kpts = []
        #start at second kpt and go to second to last segment
        nsegments = len(self.kpt_path) - 1
        for i in range(nsegments-1):
            
            #get number of points on path. this counts the first point
            try:
                i_npt = npoints[i]
            except TypeError:
                i_npt = npoints

            #this is the vector connecting the two endpoint kpts of a segment
            kdiff = self.kpt_path[i+1] - self.kpt_path[i]

            #make a vector of evenly spaced intervals, one longer than needed
            #because we chop off the last entry.
            for j in np.linspace(0,1,i_npt+1)[0:-1]:
                k = self.kpt_path[i] + j*kdiff
                #shift by small random amount to break symmetry and
                #prevent time-inversion reduction
                krand = (1. + np.random.random(3))/1.e4
                
                k += krand
                kpts.append(k)

        #now fill in the last segment, and end on the last point
        try:
            i_npt = npoints[-1]
        except TypeError:
            i_npt = npoints

        kdiff = self.kpt_path[-1] - self.kpt_path[-2]
        for j in np.linspace(0,1,i_npt+1)[1:]:
            k = self.kpt_path[-2] + j*kdiff
            #shift by small random amount to break symmetry and
            #prevent time-inversion reduction
            krand = (1. + np.random.random(3))/1.e4
            k += krand
            kpts.append(k)

        #these are now the points needed for the Harris calculation.
        self.kpts = kpts

        self.dos = DOS(self.calc)
        self.dos_energies = self.dos.get_energies()
        self.dos_dos = self.dos.get_dos()

        #try to avoid rerunning the calculation if it is already done!
        if os.path.exists(outnc):
            self.calc = Jacapo(outnc)
        else:
            print 'calculation of harris required'
            self.calc.set_nc(outnc)
            #self.calc.debug=10

            #save some time by not calculating stress
            self.calc.set_stress(False)
                        
            #this seems to be necessary sometimes
            self.calc.delete_ncattdimvar(outnc,
                                         ncdims=['number_plane_waves'])

            #this has to come after removing number_of_planewaves
            self.calc.set_kpts(self.kpts)
            
            #freeze charge density
            self.calc.set_charge_mixing(updatecharge='No')
            #and, run calculation
            self.calc.calculate()

        

    def plot(self):
        '''
        Make an interactive band-structure plot.

        clicking on a band will make it thicker and print which band was selected.
        '''

        kpoints = self.calc.get_ibz_kpoints()
        
        eigenvalues = self.calc.get_all_eigenvalues() - self.ef
        #eigenvalues = np.array([self.calc.get_eigenvalues(kpt=i)-self.ef
        #                        for i in range(len(kpoints))])
        
        self.handles = [] #used to get band indexes from plot

        fig = plt.figure()
        #plot DOS in figure
        ax = fig.add_subplot(122)
        ax.plot(self.dos_dos,self.dos_energies)
        plt.title('self-consistent Total DOS')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_ylim([-20,20])
        
        ax = fig.add_subplot(121)
        ax.set_title('Band structure')

        def onpick(event):
            'make picked line bolder, set oldline back to regular thickness'
            self.lastartist.set_linewidth(1)
            self.lastartist = thisline = event.artist
            thisline.set_linewidth(5)
            plt.draw() #needed to update linewidth
            
            print 'Band %i selected' % self.handles.index(thisline)
            #you could insert code here to plot wavefunction, etc...
            
        fig.canvas.mpl_connect('pick_event',onpick)

        #we use indices for x. the tick labels are not shown and the distance
        #appears unimportant
        xdata = range(len(eigenvalues))

        nkpts, nbands = eigenvalues.shape
        for i in range(nbands):         
            #eigenvalues has shape(nkpts,nbands)
            #note the comma after line_handle
            line_handle, = ax.plot(xdata,eigenvalues[:,i],'.-',ms=1,picker=2)
            self.handles.append(line_handle)

        self.lastartist = self.handles[-1]
            
        #plot Fermi level
        ax.plot([0,len(self.kpts)],[0,0],'k--',label='$E_f$')
        
        plt.xlabel('|k|')
        plt.ylabel('$E-E_f$ (eV)')

        #set xtick locations and labels
        xtick_locs = np.zeros(len(self.kpt_path))
        try:
            #this means the npoints is a list
            i_npt = self.npoints[0]
            for j,npt in enumerate(1,self.npoints):
                xtick_locs[j] = xtick_locs[j-1] + npt
        except TypeError:
            #npoints is a single number
            for j in range(1,len(self.labels)):
                xtick_locs[j] = xtick_locs[j-1] + self.npoints

        #the last location is off by one, so we fix it.
        xtick_locs[-1] -= 1

        ax.set_xlim([xtick_locs[0],xtick_locs[-1]])        
        ax.set_xticks(xtick_locs)
        ax.set_xticklabels(self.labels)
        
        #this seems reasonable to avoid very deep energy states and high energy states
        ax.set_ylim([-20,20])
        
        plt.show()

        return fig

 
                              

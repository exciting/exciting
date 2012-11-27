'''
python module for ASE2-free and Numeric-free dacapo

U{John Kitchin<mailto:jkitchin@andrew.cmu.edu>} December 25, 2008

This module supports numpy directly.

* ScientificPython2.8 is required

 - this is the first version to use numpy by default.

see https://wiki.fysik.dtu.dk/stuff/nc/ for dacapo netcdf variable
documentation
'''

__docformat__ = 'restructuredtext'

import sys
import exceptions, glob, os, pickle, string
from Scientific.IO.NetCDF import NetCDFFile as netCDF
import numpy as np
import subprocess as sp

import validate
import changed 

try:
    from uuid import uuid1
except ImportError: #probably an old python before 2.5
    import random, time
    def uuid1():
        t = time.asctime()
        host = os.environ['HOSTNAME']
        random.seed(host + str(t))
        s = host + '-' + t + '-'+str(random.random())
        return s.replace(' ','-')
    
import logging
log = logging.getLogger('Jacapo')

handler = logging.StreamHandler()
if sys.version_info < (2,5): # no funcName in python 2.4
    formatstring = ('%(levelname)-10s '
                    'lineno: %(lineno)-4d %(message)s')
else:
    formatstring = ('%(levelname)-10s function: %(funcName)s '
                    'lineno: %(lineno)-4d %(message)s')
formatter = logging.Formatter(formatstring)
handler.setFormatter(formatter)
log.addHandler(handler)

from ase.calculators.jacapo.validate import get_dacapopath

class DacapoRunning(exceptions.Exception):
    """Raised when ncfile.status = 'running'"""
    pass

class DacapoAborted(exceptions.Exception):
    """Raised when ncfile.status = 'aborted'"""
    pass

class DacapoInput(exceptions.Exception):
    ''' raised for bad input variables'''
    pass

class DacapoAbnormalTermination(exceptions.Exception):
    """Raised when text file does not end correctly"""
    pass

class DacapoDryrun(exceptions.Exception):
    """Raised when text file does not end correctly"""
    pass



def read(ncfile):
    '''return atoms and calculator from ncfile

    >>> atoms, calc = read('co.nc')
    '''
    calc = Jacapo(ncfile)
    atoms = calc.get_atoms() #this returns a copy
    return (atoms, calc)

class Jacapo:
    '''
    Python interface to the Fortran DACAPO code
    '''
    
    __name__ = 'Jacapo'
    __version__ = '0.4'
    
    #dictionary of valid input variables and default settings
    default_input = {'atoms':None,
                     'pw':350,
                     'dw':350,
                     'xc':'PW91',
                     'nbands':None,
                     'ft':0.1,
                     'kpts':(1,1,1),
                     'spinpol':False,
                     'fixmagmom':None,
                     'symmetry':False,
                     'calculate_stress':False,
                     'dipole':{'status':False,
                               'mixpar':0.2,
                               'initval':0.0,
                               'adddipfield':0.0,
                               'position':None},                         
                     'status':'new',
                     'pseudopotentials':None,
                     'extracharge':None,
                     'extpot':None,
                     'fftgrid':None,
                     'ascii_debug':'Off',
                     'ncoutput':{'wf':'Yes',
                                 'cd':'Yes',
                                 'efp':'Yes',
                                 'esp':'Yes'},
                     'ados':None,
                     'decoupling':None,
                     'external_dipole':None,
                     'convergence':{'energy':0.00001,
                                    'density':0.0001,
                                    'occupation':0.001,
                                    'maxsteps':None,
                                    'maxtime':None},
                     'charge_mixing':{'method':'Pulay',
                                      'mixinghistory':10,
                                      'mixingcoeff':0.1,
                                      'precondition':'No',
                                      'updatecharge':'Yes'},
                     'electronic_minimization':{'method':'eigsolve',
                                                'diagsperband':2},
                     'occupationstatistics':'FermiDirac',
                     'fftgrid':{'soft':None,
                                'hard':None},
                     'mdos':None,
                     'psp':None
                   }
    
    def __init__(self,
                 nc='out.nc',
                 outnc=None,
                 debug=logging.WARN,
                 stay_alive=False,
                 **kwargs):
        '''
        Initialize the Jacapo calculator

        :Parameters:

          nc : string
           output netcdf file, or input file if nc already exists

          outnc : string
           output file. by default equal to nc

          debug : integer
           logging debug level.
           
        Valid kwargs:

          atoms : ASE.Atoms instance
            atoms is an ase.Atoms object that will be attached
            to this calculator.
        
          pw : integer
            sets planewave cutoff

          dw : integer
            sets density cutoff

          kpts : iterable
            set chadi-cohen, monkhorst-pack kpt grid,
            e.g. kpts = (2,2,1) or explicit list of kpts
            
          spinpol : Boolean
            sets whether spin-polarization is used or not.

          fixmagmom : float
            set the magnetic moment of the unit cell. only used
            in spin polarize calculations

          ft : float
            set the Fermi temperature used in occupation smearing

          xc : string
            set the exchange-correlation functional.
            one of ['PZ','VWN','PW91','PBE','RPBE','revPBE'],

          dipole
            boolean
            turn the dipole correction on (True) or off (False)

            or:
            dictionary of parameters to fine-tune behavior
            {'status':False,
            'mixpar':0.2,
            'initval':0.0,
            'adddipfield':0.0,
            'position':None}
            
          nbands : integer
            set the number of bands

          symmetry : Boolean
            Turn symmetry reduction on (True) or off (False)

          stress : Boolean
            Turn stress calculation on (True) or off (False)

          debug : level for logging
            could be something like
            logging.DEBUG or an integer 0-50. The higher the integer,
            the less information you see set debug level (0 = off, 10 =
            extreme)
        
        Modification of the nc file only occurs at calculate time if needed

        >>> calc = Jacapo('CO.nc')
        
        reads the calculator from CO.nc if it exists or
        minimally initializes CO.nc with dimensions if it does not exist. 

        >>> calc = Jacapo('CO.nc', pw=300)

        reads the calculator from CO.nc or initializes it if
        it does not exist and changes the planewave cutoff energy to
        300eV

         >>> atoms = Jacapo.read_atoms('CO.nc')

        returns the atoms in the netcdffile CO.nc, with the calculator
        attached to it.

        >>> atoms, calc = read('CO.nc')
        
        '''
        self.debug = debug
        log.setLevel(debug)

        self.pars = Jacapo.default_input.copy()
        self.pars_uptodate = {}

        log.debug(self.pars)
        
        for key in self.pars:
            self.pars_uptodate[key] = False
            
        self.kwargs = kwargs
        self.set_psp_database()

        self.set_nc(nc)
        #assume not ready at init, rely on code later to change this 
        self.ready = False  

        # need to set a default value for stay_alive
        self.stay_alive = stay_alive

        # for correct updating, we need to set the correct frame number
        # before setting atoms or calculator
        self._set_frame_number()
        
        if os.path.exists(nc):

            self.atoms = self.read_only_atoms(nc)

            #if atoms object is passed to
            #__init__ we assume the user wants the atoms object
            # updated to the current state in the file.
            if 'atoms' in kwargs:
                log.debug('Updating the atoms in kwargs')

                atoms = kwargs['atoms']
                atoms.set_cell(self.atoms.get_cell())
                atoms.set_positions(self.atoms.get_positions())
                atoms.calc = self
                
            #update the parameter list from the ncfile
            self.update_input_parameters()
    
            self.ready = True
                   
        #change output file if needed
        if outnc:
            self.set_nc(outnc)          
            
        if len(kwargs) > 0:

            if 'stress' in kwargs:
                raise DacapoInput, '''\
                stress keyword is deprecated.
                you must use calculate_stress instead'''
            
            #make sure to set calculator on atoms if it was in kwargs
            #and do this first, since some parameters need info from atoms
            if 'atoms' in kwargs:
                #we need to set_atoms here so the atoms are written to
                #the ncfile
                self.set_atoms(kwargs['atoms'])
                kwargs['atoms'].calc = self
                del kwargs['atoms'] #so we don't call it in the next
                                    #line. we don't want to do that
                                    #because it will update the _frame
                                    #counter, and that should not be
                                    #done here.
                                
            self.set(**kwargs) #if nothing changes, nothing will be done
                    
    def set(self, **kwargs):
        '''set a parameter

        parameter is stored in dictionary that is processed later if a
        calculation is need.
        '''

        if 'DACAPO_NOSET' in os.environ:
            #it is probably a bug that this is detected so we raise an exception
            raise Exception, 'DACAPO_NOSET detected, nothing is being set'
            
        
        for key in kwargs:
            if key not in self.default_input:
                raise DacapoInput, '%s is not valid input' % key

            if kwargs[key] is None:
                continue
            
            #now check for valid input
            validf = 'validate.valid_%s' % key
            valid = eval('%s(kwargs[key])' % validf)
            if not valid:
                s = 'Warning invalid input detected for key "%s" %s'
                log.warn(s % (key,
                              kwargs[key]))
                raise DacapoInput, s % (key, kwargs[key])
                                       
            #now see if key has changed
            if key in self.pars:
                changef = 'changed.%s_changed' % key
                if os.path.exists(self.get_nc()):
                    notchanged = not eval('%s(self,kwargs[key])' % changef)
                else:
                    notchanged = False
                log.debug('%s notchanged = %s' % (key, notchanged))

                if notchanged:
                    continue

            log.debug('setting: %s. self.ready = False ' % key)

            # psp's are not stored in self.pars, everything else is
            if key == 'psp':
                self.psp[kwargs[key]['sym']] = kwargs[key]['psp']
            else:
                self.pars[key] = kwargs[key]
            self.pars_uptodate[key] = False
            self.ready = False
            log.debug('exiting set function')

    def write_input(self):
        '''write out input parameters as needed

        you must define a self._set_keyword function that does all the
        actual writing.
        '''
        
        log.debug('Writing input variables out')
        log.debug(self.pars)
        
        if 'DACAPO_READONLY' in os.environ:
            raise Exception, 'DACAPO_READONLY set and you tried to write!'
        
        if self.ready:
            log.debug('self.ready = %s' % self.ready)
            log.debug('detected everything is ready, not writing input out')
            return

        # Only write out changed parameters. this function does not do
        # the writing, that is done for each variable in private
        # functions.
        for key in self.pars:
            if self.pars_uptodate[key] is False:
                setf = 'set_%s' % key

                if self.pars[key] is None:
                    continue
                
                log.debug('trying to call: %s' % setf)
                log.debug('self.%s(self.pars[key])' % setf)
                log.debug('key = %s' % str(self.pars[key]))

                if isinstance(self.pars[key], dict):
                    eval('self.%s(**self.pars[key])' % setf)
                else:
                    eval('self.%s(self.pars[key])' % setf)
                
                self.pars_uptodate[key] = True #update the changed flag

                log.debug('wrote %s: %s' % (key, str(self.pars[key])))

        #set Jacapo version
        ncf = netCDF(self.get_nc(), 'a')
        ncf.Jacapo_version = Jacapo.__version__
        ncf.sync()
        ncf.close()

    def update_input_parameters(self):
        '''read in all the input parameters from the netcdfile'''

        log.debug('Updating parameters')
        
        for key in self.default_input:
            getf = 'self.get_%s()' % key
            log.debug('getting key: %s' % key)
            self.pars[key] = eval(getf)
            self.pars_uptodate[key] = True
        return self.pars

    def write(self, new=False):
        '''write out everything to the ncfile : self.get_nc()

        new determines whether to delete any existing ncfile, and rewrite it.
        '''
        nc = self.get_nc()

        if new:
            if os.path.exists(nc):
                os.unlink(nc)
            self.ready = False
            for key in self.pars_uptodate:
                self.pars_uptodate[key] = False

        if not os.path.exists(nc):
            self.initnc()

        self.write_input()
        self.write_nc()
        
    
    def initnc(self, ncfile=None):
        '''create an ncfile with minimal dimensions in it

        this makes sure the dimensions needed for other set functions
        exist when needed.'''
        
        if ncfile is None:
            ncfile = self.get_nc()
        else:
            self.set_nc(ncfile)
            
        log.debug('initializing %s' % ncfile)
                
        base = os.path.split(ncfile)[0]
        if base is not '' and not os.path.isdir(base):
            os.makedirs(base)
        
        ncf = netCDF(ncfile, 'w')
        #first, we define some dimensions we always need
        #unlimited
        ncf.createDimension('number_ionic_steps', None)
        ncf.createDimension('dim1', 1)
        ncf.createDimension('dim2', 2)
        ncf.createDimension('dim3', 3)
        ncf.createDimension('dim4', 4)
        ncf.createDimension('dim5', 5)
        ncf.createDimension('dim6', 6)
        ncf.createDimension('dim7', 7)
        ncf.createDimension('dim20', 20) #for longer strings
        ncf.status  = 'new'
        ncf.history = 'Dacapo'
        ncf.uuid = str(uuid1())
        ncf.Jacapo_version = Jacapo.__version__
        ncf.close()
        
        self.ready = False
        self._frame = 0
        
    def __del__(self):
        '''If calculator is deleted try to stop dacapo program
        '''
        
        if hasattr(self, '_dacapo'):
            if self._dacapo.poll()==None:
                self.execute_external_dynamics(stopprogram=True)
        #and clean up after Dacapo
        if os.path.exists('stop'):
            os.remove('stop')
        #remove slave files
        txt = self.get_txt()
        if txt is not None:
            slv = txt + '.slave*'
            for slvf in glob.glob(slv):
                os.remove(slvf)
        
    def __str__(self):
        '''
        pretty-print the calculator and atoms.

        we read everything directly from the ncfile to prevent
        triggering any calculations
        '''
        s = []
        if self.nc is None:
            return 'No netcdf file attached to this calculator'
        if not os.path.exists(self.nc):
            return 'ncfile (%s) does not exist yet' % self.nc
        
        nc = netCDF(self.nc, 'r')
        s.append('  ---------------------------------')
        s.append('  Dacapo calculation from %s' % self.nc)
        if hasattr(nc, 'uuid'):
            s.append('  uuid = %s' % nc.uuid)
        if hasattr(nc, 'status'):
            s.append('  status = %s' % nc.status)
        if hasattr(nc, 'version'):
            s.append('  version = %s' % nc.version)
        if hasattr(nc, 'Jacapo_version'):
            s.append('  Jacapo version = %s' % nc.Jacapo_version[0])
        
        energy = nc.variables.get('TotalEnergy', None)
        
        if energy and energy[:][-1] < 1E36: # missing values get
                                              # returned at 9.3E36
            s.append('  Energy = %1.6f eV' % energy[:][-1])
        else:
            s.append('  Energy = None')
            
        s.append('')
        
        atoms = self.get_atoms()
        
        if atoms is None:
            s.append('  no atoms defined')
        else:
            uc = atoms.get_cell()
            #a, b, c = uc
            s.append("  Unit Cell vectors (angstroms)")
            s.append("         x        y       z   length")

            for i, v in enumerate(uc):
                L = (np.sum(v**2))**0.5 #vector length
                s.append("  a%i [% 1.4f % 1.4f % 1.4f] %1.2f" % (i,
                                                                 v[0],
                                                                 v[1],
                                                                 v[2],
                                                                 L))
                                                                 
            stress = nc.variables.get('TotalStress', None)
            if stress is not None:
                stress = np.take(stress[:].ravel(), [0, 4, 8, 5, 2, 1])
                s.append('  Stress: xx,   yy,    zz,    yz,    xz,    xy')
                s1 = '       % 1.3f % 1.3f % 1.3f % 1.3f % 1.3f % 1.3f'
                s.append(s1 % tuple(stress))
            else:
                s.append('  No stress calculated.')
            s.append('  Volume = %1.2f A^3' % atoms.get_volume())
            s.append('')

            z = "  Atom,  sym, position (in x,y,z),     tag, rmsForce and psp"
            s.append(z)

            #this is just the ncvariable
            forces = nc.variables.get('DynamicAtomForces', None)
           
            for i, atom in enumerate(atoms):
                sym = atom.symbol
                pos = atom.position
                tag = atom.tag
                if forces is not None and (forces[:][-1][i] < 1E36).all():
                    f = forces[:][-1][i]
                    # Lars Grabow: this seems to work right for some
                    # reason, but I would expect this to be the right
                    # index order f=forces[-1][i][:]
                    # frame,atom,direction
                    rmsforce = (np.sum(f**2))**0.5
                else:
                    rmsforce = None
                    
                st = "  %2i   %3.12s  " % (i, sym)
                st += "[% 7.3f%7.3f% 7.3f] " % tuple(pos)
                st += " %2s  " % tag
                if rmsforce is not None:
                    st += " %4.3f " % rmsforce
                else:
                    st += ' None '
                st += " %s"  % (self.get_psp(sym))        
                s.append(st)
                
        s.append('')
        s.append('  Details:')
        xc = self.get_xc()
        if xc is not None:
            s.append('  XCfunctional        = %s' % self.get_xc())
        else:
            s.append('  XCfunctional        = Not defined')

        pw = self.get_pw()
        if pw is None:
            pw = 'default (350eV)'
            
        s.append('  Planewavecutoff     = %s eV' % pw)
        dw = self.get_dw()
        if dw:
            s.append('  Densitywavecutoff   = %i eV' % int(self.get_dw()))
        else:
            s.append('  Densitywavecutoff   = None')
        ft = self.get_ft()
        if ft is not None:
            s.append('  FermiTemperature    = %f kT' % ft)
        else:
            s.append('  FermiTemperature    = not defined')
        try:
            nelectrons = self.get_valence()
        except:
            nelectrons = None
        if nelectrons is not None:
            s.append('  Number of electrons = %1.1f'  % nelectrons)
        else:
            s.append('  Number of electrons = N/A')
        s.append('  Number of bands     = %s'  % self.get_nbands())
        s.append('  Kpoint grid         = %s' % str(self.get_kpts_type()))
        s.append('  Spin-polarized      = %s' % self.get_spin_polarized())
#        if self.get_spin_polarized():
#           s.append('    Unit cell magnetic moment = %1.2f bohr-magnetons' % \
#                  self.get_magnetic_moment())
        s.append('  Dipole correction   = %s' % self.get_dipole())
        s.append('  Symmetry            = %s' % self.get_symmetry())
        s.append('  Constraints         = %s' % str(atoms._get_constraints()))
        s.append('  ---------------------------------')
        nc.close()
        return string.join(s, '\n')            

    #todo figure out other xc psp databases
    def set_psp_database(self, xc=None):
        '''
        get the xc-dependent psp database

        :Parameters:

         xc : string
           one of 'PW91', 'PBE', 'revPBE', 'RPBE', 'PZ'

        
        not all the databases are complete, and that means
        some psp do not exist.

        note: this function is not supported fully. only pw91 is
        imported now. Changing the xc at this point results in loading
        a nearly empty database, and I have not thought about how to
        resolve that
        '''
        
        if xc == 'PW91' or xc is None:
            from pw91_psp import defaultpseudopotentials
        else:
            log.warn('PW91 pseudopotentials are being used!')
            #todo build other xc psp databases
            from pw91_psp import defaultpseudopotentials

        self.psp = defaultpseudopotentials

    def _set_frame_number(self, frame=None):
        '''set framenumber in the netcdf file

        this is equal to the number of ionic steps dimension'''
        
        if frame is None:
            if os.path.exists(self.nc):
                nc = netCDF(self.nc, 'r')
                # nc.dimensions['number_ionic_steps'] is None
                if 'TotalEnergy' in nc.variables:
                    number_ionic_steps = nc.variables['TotalEnergy'].shape[0]
                else:
                    number_ionic_steps = nc.variables['DynamicAtomPositions'].shape[0]
                    
                frame = number_ionic_steps - 1
                nc.close()
            else:
                if hasattr(self,'atoms'):
                    frame = 1
                else:
                    #when atoms are set, the frame will be incremented
                    frame = 0
    
##            if 'TotalEnergy' in nc.variables:
##                frame = nc.variables['TotalEnergy'].shape[0]
##                # make sure the last energy is reasonable. Sometime
##                # the field is empty if the calculation ran out of
##                # walltime for example. Empty values get returned as
##                # 9.6E36.  Dacapos energies should always be negative,
##                # so if the energy is > 1E36, there is definitely
##                # something wrong and a restart is required.
##                if nc.variables.get('TotalEnergy', None)[-1] > 1E36:
##                    log.warn("Total energy > 1E36. NC file is incomplete. \
##                    calc.restart may be required")
##                    #self.restart()
            
        log.info("Current frame number is: %i" % (frame - 1))
        self._frame = frame - 1  #netCDF starts counting with 1

    def _increment_frame(self):
        'increment the framenumber'
        
        log.debug('incrementing frame')
        self._frame += 1

    def set_pw(self, pw):
        '''set the planewave cutoff.

        :Parameters:

         pw : integer
           the planewave cutoff in eV
           
        this function checks to make sure the density wave cutoff is
        greater than or equal to the planewave cutoff.'''
        
        nc = netCDF(self.nc, 'a')
        if 'PlaneWaveCutoff' in nc.variables:
            vpw = nc.variables['PlaneWaveCutoff']
            vpw.assignValue(pw)
        else:
            vpw = nc.createVariable('PlaneWaveCutoff', 'd', ('dim1',))
            vpw.assignValue(pw)

        if 'Density_WaveCutoff' in nc.variables:
            vdw = nc.variables['Density_WaveCutoff']
            dw = vdw.getValue()
            if pw > dw:
                vdw.assignValue(pw) #make them equal
        else:
            vdw = nc.createVariable('Density_WaveCutoff', 'd', ('dim1',))
            vdw.assignValue(pw) 
        nc.close()
        self.restart() #nc dimension change for number_plane_Wave dimension
        self.set_status('new')
        self.ready = False

    def set_dw(self, dw):
        '''set the density wave cutoff energy.

        :Parameters:

          dw : integer
            the density wave cutoff

        The function checks to make sure it is not less than the
        planewave cutoff.

        Density_WaveCutoff describes the kinetic energy neccesary to
        represent a wavefunction associated with the total density,
        i.e. G-vectors for which $\vert G\vert^2$ $<$
        4*Density_WaveCutoff will be used to describe the total
        density (including augmentation charge and partial core
        density). If Density_WaveCutoff is equal to PlaneWaveCutoff
        this implies that the total density is as soft as the
        wavefunctions described by the kinetic energy cutoff
        PlaneWaveCutoff. If a value of Density_WaveCutoff is specified
        (must be larger than or equal to PlaneWaveCutoff) the program
        will run using two grids, one for representing the
        wavefunction density (softgrid_dim) and one representing the
        total density (hardgrid_dim). If the density can be
        reprensented on the same grid as the wavefunction density
        Density_WaveCutoff can be chosen equal to PlaneWaveCutoff
        (default).
        '''

        pw = self.get_pw()
        if pw > dw:
            log.warn('Planewave cutoff %i is greater \
than density cutoff %i' % (pw, dw))
        
        ncf = netCDF(self.nc, 'a')
        if 'Density_WaveCutoff' in ncf.variables:
            vdw = ncf.variables['Density_WaveCutoff']
            vdw.assignValue(dw)
        else:
            vdw = ncf.createVariable('Density_WaveCutoff', 'i', ('dim1',))
            vdw.assignValue(dw)
        ncf.close()
        self.restart() #nc dimension change
        self.set_status('new')
        self.ready = False

    def set_xc(self, xc):
        '''Set the self-consistent exchange-correlation functional

        :Parameters:

         xc : string
           Must be one of 'PZ', 'VWN', 'PW91', 'PBE', 'revPBE', 'RPBE'

        Selects which density functional to use for
        exchange-correlation when performing electronic minimization
        (the electronic energy is minimized with respect to this
        selected functional) Notice that the electronic energy is also
        evaluated non-selfconsistently by DACAPO for other
        exchange-correlation functionals Recognized options :

        * "PZ" (Perdew Zunger LDA-parametrization)
        * "VWN" (Vosko Wilk Nusair LDA-parametrization)
        * "PW91" (Perdew Wang 91 GGA-parametrization)
        * "PBE" (Perdew Burke Ernzerhof GGA-parametrization)
        * "revPBE" (revised PBE/1 GGA-parametrization)
        * "RPBE" (revised PBE/2 GGA-parametrization)

        option "PZ" is not allowed for spin polarized
        calculation; use "VWN" instead.
        '''
        nc = netCDF(self.nc, 'a')
        v = 'ExcFunctional'
        if v in nc.variables:
            nc.variables[v][:] = np.array('%7s' % xc, 'c') 
        else:
            vxc = nc.createVariable('ExcFunctional', 'c', ('dim7',))
            vxc[:] = np.array('%7s' % xc, 'c')
        nc.close()
        self.set_status('new')
        self.ready = False    

    def set_nbands(self, nbands=None):
        '''Set the number of bands. a few unoccupied bands are
        recommended.

        :Parameters:

          nbands : integer
            the number of bands.
            
        if nbands = None the function returns with nothing done. At
        calculate time, if there are still no bands, they will be set
        by:

        the number of bands is calculated as
        $nbands=nvalence*0.65 + 4$
        '''
        if nbands is None:
            return

        self.delete_ncattdimvar(self.nc,
                                ncdims=['number_of_bands'],
                                ncvars=[])
                       
        nc = netCDF(self.nc, 'a')
        v = 'ElectronicBands'
        if v in nc.variables:
            vnb = nc.variables[v]
        else:
            vnb = nc.createVariable('ElectronicBands', 'c', ('dim1',))
            
        vnb.NumberOfBands = nbands
        nc.sync()
        nc.close()
        self.set_status('new')
        self.ready = False

    def set_kpts(self, kpts):
        '''
        set the kpt grid.

        Parameters:

        kpts: (n1,n2,n3) or [k1,k2,k3,...] or one of these
        chadi-cohen sets:
         
        * cc6_1x1
        * cc12_2x3
        * cc18_sq3xsq3
        * cc18_1x1
        * cc54_sq3xsq3
        * cc54_1x1
        * cc162_sq3xsq3
        * cc162_1x1
        
        (n1,n2,n3) creates an n1 x n2 x n3 monkhorst-pack grid,
        [k1,k2,k3,...] creates a kpt-grid based on the kpoints
        defined in k1,k2,k3,...
        
        There is also a possibility to have Dacapo (fortran) create
        the Kpoints in chadi-cohen or monkhorst-pack form. To do this
        you need to set the KpointSetup.gridtype attribute, and
        KpointSetup.

        KpointSetup = [3,0,0]
        KpointSetup.gridtype = 'ChadiCohen'

        KpointSetup(1)  Chadi-Cohen k-point set
        1       6 k-points 1x1
        2       18-kpoints sqrt(3)*sqrt(3)
        3       18-kpoints 1x1
        4       54-kpoints sqrt(3)*sqrt(3)
        5       54-kpoints 1x1
        6       162-kpoints 1x1
        7       12-kpoints 2x3
        8       162-kpoints 3xsqrt 3

        or
        KpointSetup = [4,4,4]
        KpointSetup.gridtype = 'MonkhorstPack'
        we do not use this functionality.
        '''

        #chadi-cohen
        if isinstance(kpts, str):
            exec('from ase.dft.kpoints import %s' % kpts)
            listofkpts = eval(kpts)
            gridtype = kpts #stored in ncfile
            #uc = self.get_atoms().get_cell()
            #listofkpts = np.dot(ccgrid,np.linalg.inv(uc.T))

        #monkhorst-pack grid
        if np.array(kpts).shape == (3,):
            from ase.dft.kpoints import monkhorst_pack
            N1, N2, N3 = kpts
            listofkpts = monkhorst_pack((N1, N2, N3))
            gridtype = 'Monkhorst-Pack %s' % str(tuple(kpts))

        #user-defined list is provided
        if len(np.array(kpts).shape) == 2:
            listofkpts = kpts
            gridtype = 'user_defined_%i_kpts' % len(kpts)  #stored in ncfile
                    
        nbzkpts = len(listofkpts)

        #we need to get dimensions stored temporarily so
        #we can delete all dimensions and variables associated with
        #kpoints before we save them back out.
        nc2 = netCDF(self.nc, 'r')
        ncdims = nc2.dimensions
        nc2.close()

        if 'number_BZ_kpoints' in ncdims:
            self.delete_ncattdimvar(self.nc,
                                    ncdims=['number_plane_waves',
                                            'number_BZ_kpoints',
                                            'number_IBZ_kpoints'])

        # now define dim and var
        nc = netCDF(self.nc, 'a')
        nc.createDimension('number_BZ_kpoints', nbzkpts)
        bv = nc.createVariable('BZKpoints', 'd', ('number_BZ_kpoints',
                                                  'dim3'))

        bv[:] = listofkpts
        bv.gridtype = gridtype
        nc.sync()
        nc.close()

        log.debug('kpts = %s' % str(self.get_kpts()))
        
        self.set_status('new')
        self.ready = False

    def atoms_are_equal(self, atoms):
        '''
        comparison of atoms to self.atoms using tolerances to account
        for float/double differences and float math.
        '''
        
        TOL = 1.0e-6 #angstroms

        a = self.atoms.arrays
        b = atoms.arrays

        #match number of atoms in cell
        lenmatch = len(atoms) == len(self.atoms)
        if lenmatch is not True:
            return False #the next two comparisons fail in this case.
        
        #match positions in cell
        posmatch = (abs(a['positions'] - b['positions']) <= TOL).all()
        #match cell
        cellmatch = (abs(self.atoms.get_cell()
                         - atoms.get_cell()) <= TOL).all()
        
        if lenmatch and posmatch and cellmatch:
            return True
        else:
            return False
                
    def set_atoms(self, atoms):
        '''attach an atoms to the calculator and update the ncfile

        :Parameters:

          atoms
            ASE.Atoms instance
          
        '''

        log.debug('setting atoms to: %s' % str(atoms))
        
        if hasattr(self, 'atoms') and self.atoms is not None:
            #return if the atoms are the same. no change needs to be made
            if self.atoms_are_equal(atoms):
                log.debug('No change to atoms in set_atoms, returning')
                return
            
            # some atoms already exist. Test if new atoms are
            # different from old atoms.
            # this is redundant
            if atoms != self.atoms:
                # the new atoms are different from the old ones. Start
                # a new frame.
                log.debug('atoms != self.atoms, incrementing')
                self._increment_frame()
                
        self.atoms = atoms.copy()
        self.ready = False
        log.debug('self.atoms = %s' % str(self.atoms))

    def set_ft(self, ft):
        '''set the Fermi temperature for occupation smearing

        :Parameters:

          ft : float
            Fermi temperature in kT (eV)

        Electronic temperature, corresponding to gaussian occupation
        statistics. Device used to stabilize the convergence towards
        the electronic ground state. Higher values stabilizes the
        convergence. Values in the range 0.1-1.0 eV are recommended,
        depending on the complexity of the Fermi surface (low values
        for d-metals and narrow gap semiconducters, higher for free
        electron-like metals).
        '''
        
        nc = netCDF(self.nc, 'a')
        v = 'ElectronicBands'
        if v in nc.variables:
            vnb = nc.variables[v]
        else:
            vnb = nc.createVariable('ElectronicBands', 'c', ('dim1',))

        vnb.OccupationStatistics_FermiTemperature = ft
        nc.sync()
        nc.close()
        self.set_status('new')
        self.ready = False

    def set_status(self, status):
        '''set the status flag in the netcdf file

        :Parameters:

          status : string
            status flag, e.g. 'new', 'finished'
        '''
        
        nc = netCDF(self.nc, 'a')
        nc.status = status
        nc.sync()
        nc.close()
        log.debug('set status to %s' % status)

    def get_spinpol(self):
        'Returns the spin polarization setting, either True or False'
        
        nc = netCDF(self.nc, 'r')
        v = 'ElectronicBands'
        if v in nc.variables:
            vnb = nc.variables[v]
            if hasattr(vnb, 'SpinPolarization'):
                spinpol = vnb.SpinPolarization
            else:
                spinpol = 1
        else:
            spinpol = 1

        nc.close()
        if spinpol == 1:
            return False
        else:
            return True

    def set_spinpol(self, spinpol=False):
        '''set Spin polarization.

        :Parameters:

         spinpol : Boolean
           set_spinpol(True)  spin-polarized.
           set_spinpol(False) no spin polarization, default

        Specify whether to perform a spin polarized or unpolarized
        calculation.
        '''
        
        nc = netCDF(self.nc, 'a')
        v = 'ElectronicBands'
        if v in nc.variables:
            vnb = nc.variables[v]
        else:
            vnb = nc.createVariable('ElectronicBands', 'c', ('dim1',))

        if spinpol is True:
            vnb.SpinPolarization = 2
        else:
            vnb.SpinPolarization = 1

        nc.sync()
        nc.close()
        self.set_status('new')
        self.ready = False

    def set_fixmagmom(self, fixmagmom=None):
        '''set a fixed magnetic moment for a spin polarized calculation

        :Parameters:

          fixmagmom : float
            the magnetic moment of the cell in Bohr magnetons
        '''
        
        if fixmagmom is None:
            return
        
        nc = netCDF(self.nc,'a')
        v = 'ElectronicBands'
        if v in nc.variables:
            vnb = nc.variables[v]
        else:
            vnb = nc.createVariable('ElectronicBands', 'c', ('dim1',))

        vnb.SpinPolarization = 2 #You must want spin-polarized
        vnb.FixedMagneticMoment = fixmagmom
        nc.sync()
        nc.close()
        self.set_status('new')
        self.ready = False

    def get_fixmagmom(self):
        'returns the value of FixedMagneticMoment'
        
        nc = netCDF(self.nc,'r')
        if 'ElectronicBands' in nc.variables:
            v = nc.variables['ElectronicBands']
            if hasattr(v,'FixedMagneticMoment'):
                fixmagmom = v.FixedMagneticMoment
            else:
                fixmagmom = None
        else:
            fixmagmom = None
        nc.close()
        return fixmagmom
            
    def set_calculate_stress(self, stress=True):
        '''Turn on stress calculation

        :Parameters:

          stress : boolean
            set_calculate_stress(True) calculates stress
            set_calculate_stress(False) do not calculate stress
        '''
        
        nc = netCDF(self.get_nc(),'a')
        vs = 'NetCDFOutputControl'
        if vs in nc.variables:
            v = nc.variables[vs]
        else:
            v = nc.createVariable('NetCDFOutputControl', 'c', ('dim1',))

        if stress is True:
            v.PrintTotalStress = 'Yes'
        else:
            v.PrintTotalStress = 'No'
        nc.sync()
        nc.close()
        self.set_status('new')
        self.ready = False

    def set_nc(self, nc='out.nc'):
        '''
        set filename for the netcdf and text output for this calculation

        :Parameters:

          nc : string
            filename for netcdf file
                
        if the ncfile attached to the calculator is changing, the old
        file will be copied to the new file if it doesn not exist so
        that all the calculator details are preserved. Otherwise, the 

        if the ncfile does not exist, it will get initialized.

        the text file will have the same basename as the ncfile, but
        with a .txt extension.
        '''
        
        #the first time this is called, there may be no self.nc defined
        if not hasattr(self, 'nc'):
            self.nc = nc    
            
        #check if the name is changing and if so, copy the old ncfile
        #to the new one.  This is necessary to ensure all the
        #calculator details are copied over. if the file already
        #exists we use the contents of the existing file
        if nc != self.nc and not os.path.exists(nc):
            log.debug('copying %s to %s' % (self.nc, nc))
            #import shutil
            #shutil.copy(self.nc,nc)
            base = os.path.split(nc)[0]
            if not os.path.isdir(base) and base is not '':
                os.makedirs(base)
            status = os.system("cp '%s' '%s'" % (self.nc, nc))
            if status != 0:
                raise Exception, 'Copying ncfile failed.'
            self.nc = nc
        
        elif os.path.exists(nc):
            self._set_frame_number()
            self.set_psp_database()
            self.atoms = self.read_only_atoms(nc)
            self.nc = nc
            self.update_input_parameters()
                
        #I always want the text file set based on the ncfile
        #and I never want to set this myself.
        base = os.path.splitext(self.nc)[0]
        self.txt = "%s.txt" % base

    def set_pseudopotentials(self, pspdict):
        '''Set all the pseudopotentials from a dictionary.

        The dictionary should have this form::

            {symbol1: path1,
             symbol2: path2}
        '''
        for key in pspdict:
            self.set_psp(sym=key,
                         psp=pspdict[key])
            
    def set_psp(self,
                sym=None,
                z=None,
                psp=None):
        '''
        set the pseudopotential file for a species or an atomic number.

        :Parameters:

         sym : string
           chemical symbol of the species

          z : integer
            the atomic number of the species

          psp : string
            filename of the pseudopotential

        
        you can only set sym or z.

        examples::
        
          set_psp('N',psp='pspfile')
          set_psp(z=6,psp='pspfile')
        '''
        log.debug(str([sym, z, psp]))
        if (sym, z, psp) == (None, None, None):
            return
        
        if (sym is None and z is not None):
            from ase.data import chemical_symbols
            sym = chemical_symbols[z]
        elif (sym is not None and z is None):
            pass
        else:
            raise Exception, 'You can only specify Z or sym!'

        if not hasattr(self, 'psp'):
            self.set_psp_database()

        #only make change if needed
        if sym not in self.psp:
            self.psp[sym] = psp
            self.ready = False
            self.set_status('new')
        elif self.psp[sym] != psp:
            self.psp[sym] = psp
            self.ready = False
            self.set_status('new')

        if not self.ready:
            #now we update the netcdf file
            ncf = netCDF(self.nc, 'a')
            vn = 'AtomProperty_%s' % sym
            if vn not in ncf.variables:
                if 'dim20' not in ncf.dimensions:
                    ncf.createDimension('dim20', 20)
                p = ncf.createVariable(vn, 'c', ('dim20',))
            else:
                p = ncf.variables[vn]

            ppath = self.get_psp(sym=sym)
            p.PspotFile = ppath
            ncf.close()

    def get_pseudopotentials(self):
        'get pseudopotentials set for atoms attached to calculator'
        
        if self.atoms is None:
            return None
        
        psp = {}
        for atom in self.atoms:
            psp[atom.symbol] = self.psp[atom.symbol]
        return {'pspdict':psp}

    def get_symmetry(self):
        '''return the type of symmetry used'''
        
        nc = netCDF(self.nc, 'r')
        if 'UseSymmetry' in nc.variables:
            sym = string.join(nc.variables['UseSymmetry'][:],'').strip()
        else:
            sym = None      
        nc.close()
        if sym in ['Off', None]:
            return False
        elif sym == 'Maximum':
            return True
        else:
            raise Exception, 'Type of symmetry not recognized: %s' % sym
       
    def set_symmetry(self, val=False):
        '''set how symmetry is used to reduce k-points

        :Parameters:

         val : Boolean
           set_sym(True) Maximum symmetry is used
           set_sym(False) No symmetry is used

        This variable controls the if and how DACAPO should attempt
        using symmetry in the calculation. Imposing symmetry generally
        speeds up the calculation and reduces numerical noise to some
        extent. Symmetry should always be applied to the maximum
        extent, when ions are not moved. When relaxing ions, however,
        the symmetry of the equilibrium state may be lower than the
        initial state. Such an equilibrium state with lower symmetry
        is missed, if symmetry is imposed. Molecular dynamics-like
        algorithms for ionic propagation will generally not break the
        symmetry of the initial state, but some algorithms, like the
        BFGS may break the symmetry of the initial state. Recognized
        options:

        "Off": No symmetry will be imposed, apart from time inversion
        symmetry in recipical space. This is utilized to reduce the
        k-point sampling set for Brillouin zone integration and has no
        influence on the ionic forces/motion.

        "Maximum": DACAPO will look for symmetry in the supplied
        atomic structure and extract the highest possible symmetry
        group. During the calculation, DACAPO will impose the found
        spatial symmetry on ionic forces and electronic structure,
        i.e. the symmetry will be conserved during the calculation.
        '''
        
        if val:
            symval = 'Maximum'
        else:
            symval = 'Off'
        
        ncf = netCDF(self.get_nc(), 'a')
        if 'UseSymmetry' not in ncf.variables:
            sym = ncf.createVariable('UseSymmetry', 'c', ('dim7',))
        else:
            sym = ncf.variables['UseSymmetry']
            
        sym[:] = np.array('%7s' % symval, 'c')
        ncf.sync()
        ncf.close()
        self.set_status('new')
        self.ready = False

    def set_extracharge(self, val):
        '''add extra charge to unit cell

        :Parameters:

          val : float
            extra electrons to add or subtract from the unit cell

        Fixed extra charge in the unit cell (i.e. deviation from
        charge neutrality). This assumes a compensating, positive
        constant backgound charge (jellium) to forge overall charge
        neutrality.
        '''
        
        nc = netCDF(self.get_nc(), 'a')
        if 'ExtraCharge' in nc.variables:
            v = nc.variables['ExtraCharge']
        else:
            v = nc.createVariable('ExtraCharge', 'd', ('dim1',))

        v.assignValue(val)
        nc.sync()
        nc.close()

    def get_extracharge(self):
        'Return the extra charge set in teh calculator'
        
        nc = netCDF(self.get_nc(), 'r')
        if 'ExtraCharge' in nc.variables:
            v = nc.variables['ExtraCharge']
            exchg = v.getValue()
        else:
            exchg = None
        nc.close()
        return exchg

    def get_extpot(self):
        'return the external potential set in teh calculator'
        
        nc = netCDF(self.get_nc(), 'r')
        if 'ExternalPotential' in nc.variables:
            v = nc.variables['ExternalPotential']
            extpot = v[:]
        else:
            extpot = None

        nc.close()
        return extpot
    
    def set_extpot(self, potgrid):
        '''add external potential of value

        see this link before using this
        https://listserv.fysik.dtu.dk/pipermail/campos/2003-August/000657.html
        
        :Parameters:

          potgrid : np.array with shape (nx,ny,nz)
            the shape must be the same as the fft soft grid
            the value of the potential to add

        
        you have to know both of the fft grid dimensions ahead of time!
        if you know what you are doing, you can set the fft_grid you want
        before hand with:
        calc.set_fftgrid((n1,n2,n3))
        '''
        
        nc = netCDF(self.get_nc(), 'a')
        if 'ExternalPotential' in nc.variables:
            v = nc.variables['ExternalPotential']
        else:
            # I assume here you have the dimensions of potgrid correct
            # and that the soft and hard grids are the same. 
            # if softgrid is defined, Dacapo requires hardgrid to be
            # defined too.
            s1, s2, s3 = potgrid.shape
            if 'softgrid_dim1' not in nc.dimensions:
                nc.createDimension('softgrid_dim1', s1)
                nc.createDimension('softgrid_dim2', s2)
                nc.createDimension('softgrid_dim3', s3)
                nc.createDimension('hardgrid_dim1', s1)
                nc.createDimension('hardgrid_dim2', s2)
                nc.createDimension('hardgrid_dim3', s3)
                
            v = nc.createVariable('ExternalPotential',
                                  'd',
                                  ('softgrid_dim1',
                                   'softgrid_dim2',
                                   'softgrid_dim3',))
        v[:] = potgrid
        nc.sync()
        nc.close()
        self.set_status('new')
        self.ready = False
        
    def set_fftgrid(self, soft=None, hard=None):
        '''
        sets the dimensions of the FFT grid to be used

        :Parameters:

          soft : (n1,n2,n3) integers
            make a n1 x n2 x n3 grid

          hard : (n1,n2,n3) integers
            make a n1 x n2 x n3 grid

        
        >>> calc.set_fftgrid(soft=[42,44,46])
        sets the soft and hard grid dimensions to 42,44,46

        >>> calc.set_fftgrid(soft=[42,44,46],hard=[80,84,88])
        sets the soft grid dimensions to 42,44,46 and the hard
        grid dimensions to 80,84,88
        
        These are the fast FFt grid numbers listed in fftdimensions.F
        
        data list_of_fft /2,  4,  6,  8, 10, 12, 14, 16, 18, 20, &
        22,24, 28, 30,32, 36, 40, 42, 44, 48, &
        56,60, 64, 66, 70, 72, 80, 84, 88, 90, &
        96,108,110,112,120,126,128,132,140,144,154, &
        160,168,176,180,192,198,200, &
        216,240,264,270,280,288,324,352,360,378,384,400,432, &
        450,480,540,576,640/
        
        otherwise you will get some errors from mis-dimensioned variables.
        
        this is usually automatically set by Dacapo.
        '''
        
        if soft is not None:
            self.delete_ncattdimvar(self.nc,
                                    ncdims=['softgrid_dim1',
                                            'softgrid_dim2',
                                            'softgrid_dim3'
                                            ],
                                    ncvars=[])
        
            
            nc = netCDF(self.get_nc(), 'a')
            nc.createDimension('softgrid_dim1', soft[0])
            nc.createDimension('softgrid_dim2', soft[1])
            nc.createDimension('softgrid_dim3', soft[2])
            nc.sync()
            nc.close()

            if hard is None:
                hard = soft

        if hard is not None:
            self.delete_ncattdimvar(self.nc,
                                    ncdims=['hardgrid_dim1',
                                            'hardgrid_dim2',
                                            'hardgrid_dim3'
                                            ],
                                    ncvars=[])
            nc = netCDF(self.get_nc(),'a')
            nc.createDimension('hardgrid_dim1', hard[0])
            nc.createDimension('hardgrid_dim2', hard[1])
            nc.createDimension('hardgrid_dim3', hard[2])
            nc.sync()
            nc.close()

        self.set_status('new')
        self.ready = False

    def get_ascii_debug(self):
        'Return the debug settings in Dacapo'
        
        nc = netCDF(self.get_nc(), 'r')
        if 'PrintDebugInfo' in nc.variables:
            v = nc.variables['PrintDebugInfo']
            debug = string.join(v[:], '')
        else:
            debug = None
        nc.close()
        return debug
    
    def set_ascii_debug(self, level):
        '''set the debug level for Dacapo

        :Parameters:
        
          level : string
            one of 'Off', 'MediumLevel', 'HighLevel'
        '''
        
        nc = netCDF(self.get_nc(), 'a')
        if 'PrintDebugInfo' in nc.variables:
            v = nc.variables['PrintDebugInfo']
        else:
            if 'dim20' not in nc.dimensions:
                nc.createDimension('dim20', 20)
            v = nc.createVariable('PrintDebugInfo', 'c', ('dim20',))

        v[:] = np.array('%20s' % level, dtype='c')
        nc.sync()
        nc.close()
        self.set_status('new')
        self.ready = False

    def get_ncoutput(self):
        'returns the control variables for the ncfile'
        
        nc = netCDF(self.get_nc(), 'r')
        if 'NetCDFOutputControl' in nc.variables:
            v = nc.variables['NetCDFOutputControl']
            ncoutput = {}
            if hasattr(v, 'PrintWaveFunction'):
                ncoutput['wf'] = v.PrintWaveFunction
            if hasattr(v, 'PrintChargeDensity'):
                ncoutput['cd'] = v.PrintChargeDensity
            if hasattr(v, 'PrintEffPotential'):
                ncoutput['efp'] = v.PrintEffPotential
            if hasattr(v, 'PrintElsPotential'):
                ncoutput['esp'] = v.PrintElsPotential
        else:
            ncoutput = None
        nc.close()
        return ncoutput
        
    def set_ncoutput(self,
                     wf=None,
                     cd=None,
                     efp=None,
                     esp=None):
        '''set the output of large variables in the netcdf output file

        :Parameters:

          wf : string
            controls output of wavefunction. values can
            be 'Yes' or 'No'

          cd : string
            controls output of charge density. values can
            be 'Yes' or 'No'

          efp : string
            controls output of effective potential. values can
            be 'Yes' or 'No'

          esp : string
            controls output of electrostatic potential. values can
            be 'Yes' or 'No'
        '''
        nc = netCDF(self.get_nc(), 'a')
        if 'NetCDFOutputControl' in nc.variables:
            v = nc.variables['NetCDFOutputControl']
        else:
            v = nc.createVariable('NetCDFOutputControl', 'c', ())

        if wf is not None:
            v.PrintWaveFunction = wf
        if cd is not None:
            v.PrintChargeDensity = cd
        if efp is not None:
            v.PrintEffPotential = efp
        if esp is not None:
            v.PrintElsPotential = esp

        nc.sync()
        nc.close()
        self.set_status('new')
        self.ready = False

    def get_ados(self, **kwargs):
        '''
        attempt at maintaining backward compatibility with get_ados
        returning data

        Now when we call calc.get_ados() it will return settings,

        and calc.get_ados(atoms=[],...) should return data

        '''
        
        if len(kwargs) != 0:
            return self.get_ados_data(**kwargs)
        
        nc = netCDF(self.get_nc(),'r')
        if 'PrintAtomProjectedDOS' in nc.variables:
            v = nc.variables['PrintAtomProjectedDOS']
            ados = {}
            if hasattr(v, 'EnergyWindow'):
                ados['energywindow'] = v.EnergyWindow
            if hasattr(v, 'EnergyWidth'):
                ados['energywidth'] = v.EnergyWidth[0]
            if hasattr(v, 'NumberEnergyPoints'):
                ados['npoints'] = v.NumberEnergyPoints[0]
            if hasattr(v, 'CutoffRadius'):
                ados['cutoff'] = v.CutoffRadius[0]
        else:
            ados = None

        nc.close()
        return ados
        
    def set_ados(self,
                 energywindow=(-15,5),
                 energywidth=0.2,
                 npoints=250,
                 cutoff=1.0):
        '''
        setup calculation of atom-projected density of states

        :Parameters:

          energywindow : (float, float)
            sets (emin,emax) in eV referenced to the Fermi level

          energywidth : float
            the gaussian used in smearing

          npoints : integer
            the number of points to sample the DOS at

          cutoff : float
            the cutoff radius in angstroms for the integration.
        '''
        
        nc = netCDF(self.get_nc(), 'a')
        if 'PrintAtomProjectedDOS' in nc.variables:
            v = nc.variables['PrintAtomProjectedDOS']
        else:
            v = nc.createVariable('PrintAtomProjectedDOS', 'c', ())

        v.EnergyWindow = energywindow
        v.EnergyWidth  = energywidth
        v.NumberEnergyPoints = npoints
        v.CutoffRadius = cutoff

        nc.sync()
        nc.close()
        self.set_status('new')
        self.ready = False

    def get_mdos(self):
        'return multicentered projected dos parameters'
        nc = netCDF(self.get_nc(),'r')

        mdos = {}
        
        if 'MultiCenterProjectedDOS' in nc.variables:
            v = nc.variables['MultiCenterProjectedDOS']
            mdos['energywindow'] = v.EnergyWindow
            mdos['energywidth'] = v.EnergyWidth
            mdos['numberenergypoints'] = v.NumberEnergyPoints
            mdos['cutoffradius'] = v.CutoffRadius
            mdos['mcenters'] = eval(v.mcenters)
                
        nc.close()

        return mdos

    def get_mdos_data(self,
                      spin=0,
                      cutoffradius='infinite'):
        '''returns data from multicentered projection


        returns (mdos, rotmat)

        the rotation matrices are retrieved from the text file. I am
        not sure what you would do with these, but there was a note
        about them in the old documentation so I put the code to
        retrieve them here. the syntax for the return value is:
        rotmat[atom#][label] returns the rotation matrix for the
        center on the atom# for label. I do not not know what the
        label refers to.
        '''
        
        if self.calculation_required():
            self.calculate()
        
        nc = netCDF(self.get_nc(),'r')
        icut = 1 #short
        if cutoffradius == "infinite":
            icut = 0
            
        #var = nc.variables['MultiCenterProjectedDOS']
        integrated = nc.variables['MultiCenterProjectedDOS_IntegratedDOS'][:]
        tz = 'MultiCenterProjectedDOS_EnergyResolvedDOS'
        energyresolved = nc.variables[tz][:]
        energygrid = nc.variables['MultiCenterProjectedDOS_EnergyGrid'][:]

        number_of_multicenters  = integrated.shape[0]
        #number_of_cutoff = integrated.shape[1]
        #number_of_spin = integrated.shape[2]

        multicenterprojections = []
        for multicenter in range(number_of_multicenters): 
            #orbitals = var[multicenter]
            energyresolveddata = energyresolved[multicenter, icut, spin, :]
            #integrateddata     = integrated[multicenter, icut, spin]
            multicenterprojections.append([energygrid, energyresolveddata])

        log.info('Found %d multicenters' % len(multicenterprojections))
        nc.close()

        #now parse the text file for the rotation matrices
        rot_mat_lines = []
        txt = self.get_txt()
        if os.path.exists(txt):
            f = open(txt,'r')
            for line in f:
                if 'MUL: Rmatrix' in line:
                    rot_mat_lines.append(line)
            f.close()

            rotmat = []
            for line in rot_mat_lines:
                fields = line.split()
                novl = int(fields[2])
                ncen = int(fields[3])
                row = [float(x) for x in fields[4:]]

                try:
                    rotmat[novl-1][ncen-1].append(row)
                except IndexError:
                    try:
                        rotmat[novl-1].append([])
                        rotmat[novl-1][ncen-1].append(row)
                    except IndexError:
                        rotmat.append([])
                        rotmat[novl-1].append([])
                    rotmat[novl-1][ncen-1].append(row)
        else:
            rotmat = None
                    
        return (multicenterprojections, rotmat)
            
    def set_mdos(self,
                 mcenters=None,
                 energywindow=(-15,5),
                 energywidth=0.2,
                 numberenergypoints=250,
                 cutoffradius=1.0):
        '''Setup multicentered projected DOS.

        mcenters
           a list of tuples containing (atom#,l,m,weight)
           (0,0,0,1.0) specifies (atom 0, l=0, m=0, weight=1.0) an s orbital
           on atom 0
           
           (1,1,1,1.0) specifies (atom 1, l=1, m=1, weight=1.0) a p orbital
           with m = +1 on atom 0

           l=0 s-orbital
           l=1 p-orbital
           l=2 d-orbital

           m in range of ( -l ... 0 ... +l )

           The direction cosines for which the spherical harmonics are
           set up are using the next different atom in the list
           (cyclic) as direction pointer, so the z-direction is chosen
           along the direction to this next atom. At the moment the
           rotation matrices is only given in the text file, you can
           use grep'MUL: Rmatrix' out_o2.txt to get this information.
           
        adapated from old MultiCenterProjectedDOS.py
        '''
        if mcenters is None:
            return
        
        nc = netCDF(self.get_nc(), 'a')
        
        _listofmcenters_ = mcenters
        
        # get number of multi centers
        ncenters = len(_listofmcenters_)
        # get max number of orbitals any center 
        max_orbitals = max(map(len, _listofmcenters_))

        mmatrix = np.zeros([ncenters, max_orbitals, 4], np.float)
        ncenter = 0
        for multicenter in _listofmcenters_: 
            norbital = 0
            for orbital in multicenter: 
                mmatrix[ncenter, norbital] = orbital 
                norbital = norbital + 1 

            # signal that this multicenter contains less than
            # max_orbital orbitals
            if len(multicenter) < max_orbitals: 
                mmatrix[ncenter, len(multicenter):max_orbitals] = (-1.0, 0,
                                                                   0, 0)

            ncenter = ncenter + 1

        nc.createDimension('max_orbitals', max_orbitals)
        nc.createDimension('number_of_multicenters', ncenters)

        if 'MultiCenterProjectedDOS' in nc.variables:
            v = nc.variables['MultiCenterProjectedDOS']
        else:
            v = nc.createVariable('MultiCenterProjectedDOS',
                                  'd',
                                  ('number_of_multicenters',
                                   'max_orbitals',
                                   'dim4'))

        v.EnergyWindow = energywindow
        v.EnergyWidth = energywidth
        v.NumberEnergyPoints = numberenergypoints
        v.CutoffRadius = cutoffradius

        #this is kind of hacky, but it is needed for get_mdos so you
        #can tell if the input is changed.
        v.mcenters = str(mcenters)

        v[:] = mmatrix

        nc.sync()
        nc.close()        

    def set_debug(self, debug):
        '''
        set debug level for python logging

        debug should be an integer from 0-100 or one of the logging
        constants like logging.DEBUG, logging.WARN, etc...

        '''
        
        self.debug = debug
        log.setLevel(debug)

    def get_debug(self):
        'Return the python logging level'
        
        return self.debug

    def get_decoupling(self):
        'return the electrostatic decoupling parameters'
        
        nc = netCDF(self.get_nc(), 'r')
        if 'Decoupling' in nc.variables:
            v = nc.variables['Decoupling']
            decoupling = {}
            if hasattr(v,'NumberOfGaussians'):
                decoupling['ngaussians'] = v.NumberOfGaussians
            if hasattr(v,'ECutoff'):
                decoupling['ecutoff'] = v.ECutoff
            if hasattr(v,'WidthOfGaussian'):
                decoupling['gausswidth'] = v.WidthOfGaussian
        else:
            decoupling = None
        nc.close()
        return decoupling
        
    def set_decoupling(self,
                       ngaussians=3,
                       ecutoff=100,
                       gausswidth=0.35):
        '''
        Decoupling activates the three dimensional electrostatic
        decoupling. Based on paper by Peter E. Bloechl: JCP 103
        page7422 (1995).

        :Parameters:

          ngaussians : int
            The number of gaussian functions per atom
            used for constructing the model charge of the system

          ecutoff : int
            The cut off energy (eV) of system charge density in
            g-space used when mapping constructing the model change of
            the system, i.e. only charge density components below
            ECutoff enters when constructing the model change.

          gausswidth : float
            The width of the Gaussians defined by
            $widthofgaussian*1.5^(n-1)$  $n$=(1 to numberofgaussians)
            
        '''

        nc = netCDF(self.get_nc(), 'a')
        if 'Decoupling' in nc.variables:
            v = nc.variables['Decoupling']
        else:
            v = nc.createVariable('Decoupling', 'c', ())

        v.NumberOfGaussians = ngaussians
        v.ECutoff = ecutoff
        v.WidthOfGaussian = gausswidth

        nc.sync()
        nc.close()
        self.set_status('new')
        self.ready = False

    def set_external_dipole(self,
                            value,
                            position=None):
        '''
        Externally imposed dipole potential. This option overwrites
        DipoleCorrection if set. 

        :Parameters:

          value : float
            units of volts

          position : float
            scaled coordinates along third unit cell direction.
            if None, the compensation dipole layer plane in the
            vacuum position farthest from any other atoms on both
            sides of the slab. Do not set to 0.0.
        '''
        
        var = 'ExternalDipolePotential'
        nc = netCDF(self.get_nc(), 'a')
        if var in nc.variables:
            v = nc.variables[var]
        else:
            v = nc.createVariable('ExternalDipolePotential', 'd', ())
        
        v.assignValue(value)
        if position is not None:
            v.DipoleLayerPosition = position

        nc.sync()
        nc.close()
        self.set_status('new')
        self.ready = False

    def get_external_dipole(self):
        'return the External dipole settings'

        var = 'ExternalDipolePotential'
        nc = netCDF(self.get_nc(),'r')
        if var in nc.variables:
            v = nc.variables[var]
            value = v.getValue()
            if hasattr(v, 'DipoleLayerPosition'):
                position = v.DipoleLayerPosition
            else:
                position = None

            ed = {'value':value, 'position':position}
        else:
            ed = None
        nc.close()
        return ed
            
    def set_dipole(self,
                   status=True,
                   mixpar=0.2,
                   initval=0.0,
                   adddipfield=0.0,
                   position=None):
        '''turn on and set dipole correction scheme

        :Parameters:

          status : Boolean
            True turns dipole on. False turns Dipole off

          mixpar : float
            Mixing Parameter for the the dipole correction field
            during the electronic minimization process. If instabilities
            occur during electronic minimization, this value may be
            decreased.

          initval : float
            initial value to start at

          adddipfield : float
            additional dipole field to add
            units : V/ang
            External additive, constant electrostatic field along
            third unit cell vector, corresponding to an external
            dipole layer. The field discontinuity follows the position
            of the dynamical dipole correction, i.e. if
            DipoleCorrection:DipoleLayerPosition is set, the field
            discontinuity is at this value, otherwise it is at the
            vacuum position farthest from any other atoms on both
            sides of the slab.

          position : float
            scaled coordinates along third unit cell direction.
            If this attribute is set, DACAPO will position the
            compensation dipole layer plane in at the provided value.
            If this attribute is not set, DACAPO will put the compensation
            dipole layer plane in the vacuum position farthest from any
            other atoms on both sides of the slab. Do not set this to
            0.0

        
        calling set_dipole() sets all default values.
            
        '''
        if status == False:
            self.delete_ncattdimvar(self.nc, ncvars=['DipoleCorrection'])
            return
        
        ncf = netCDF(self.get_nc(), 'a')
        if 'DipoleCorrection' not in ncf.variables:
            dip = ncf.createVariable('DipoleCorrection', 'c', ())
        else:
            dip = ncf.variables['DipoleCorrection']
        dip.MixingParameter = mixpar
        dip.InitialValue = initval
        dip.AdditiveDipoleField = adddipfield

        if position is not None:
            dip.DipoleLayerPosition = position
            
        ncf.sync()
        ncf.close()
        self.set_status('new')
        self.ready = False
       
    def set_stay_alive(self, value):
        'set the stay alive setting'
        
        self.delete_ncattdimvar(self.nc,
                                ncvars=['Dynamics'])

        if (hasattr(self,'parent') or hasattr(self,'children')) and value == True:
            log.debug("This is a parent/child calculator and stay_alive must be false.")
            value = False

        if value in [True, False]:
            self.stay_alive = value
            #self._dacapo_is_running = False
        else:
            log.debug("stay_alive must be boolean. Value was not changed.")

    def get_stay_alive(self):
        'return the stay alive settings'
        
        return self.stay_alive

    def get_fftgrid(self):
        'return soft and hard fft grids'
        
        nc = netCDF(self.nc, 'r')
        soft = []
        hard = []
        for d in [1, 2, 3]:
            sd = 'softgrid_dim%i' % d
            hd = 'hardgrid_dim%i' % d
            if sd in nc.dimensions:
                soft.append(nc.dimensions[sd])
                hard.append(nc.dimensions[hd])
        nc.close()
        if soft == []:
            soft = None
        if hard == []:
            hard = None
        return ({'soft':soft,
                 'hard':hard})

    def get_kpts_type(self):
        'return the kpt grid type'
        
        nc = netCDF(self.nc, 'r')

        if 'BZKpoints' in nc.variables:
            bv = nc.variables['BZKpoints']
            if hasattr(bv, 'gridtype'):
                kpts_type = bv.gridtype #string saved in jacapo
            else:
                #no grid attribute, this ncfile was created pre-jacapo
                kpts_type = '%i kpts' % len(bv[:])
        else:
            kpts_type = 'BZKpoints not defined. [[0,0,0]] used by default.'

        nc.close()
        return kpts_type
    
    def get_kpts(self):
        'return the BZ kpts'
        nc = netCDF(self.nc, 'r')

        if 'BZKpoints' in nc.variables:
            bv = nc.variables['BZKpoints']
            kpts = bv[:]
        else:
            kpts = np.array(([0, 0, 0])) #default Gamma point used in Dacapo when
                             #BZKpoints not defined

        nc.close()
        return kpts
        
    def get_nbands(self):
        'return the number of bands used in the calculation'
        nc = netCDF(self.nc, 'r')

        if 'ElectronicBands' in nc.variables:
            v = nc.variables['ElectronicBands']
            if hasattr(v, 'NumberOfBands'):
                nbands = int(v.NumberOfBands[0])
            else:
                nbands = None
        else:
            nbands = None
            
        nc.close()
        return nbands
    
    def get_ft(self):
        'return the FermiTemperature used in the calculation'
        nc = netCDF(self.nc, 'r')

        if 'ElectronicBands' in nc.variables:
            v = nc.variables['ElectronicBands']
            if hasattr(v, 'OccupationStatistics_FermiTemperature'):
                ft = v.OccupationStatistics_FermiTemperature
            else:
                ft = None
        else:
            ft = None
        nc.close()
        return ft
    
    def get_dipole(self):
        'return dictionary of parameters if the DipoleCorrection was used'
        
        nc = netCDF(self.get_nc(), 'r')
        pars = {}
        if 'DipoleCorrection' in nc.variables:
            v = nc.variables['DipoleCorrection']
            pars['status'] = True
            if hasattr(v, 'MixingParameter'):
                pars['mixpar'] = v.MixingParameter
            if hasattr(v, 'InitialValue'):
                pars['initval'] = v.InitialValue
            if hasattr(v, 'AdditiveDipoleField'):
                pars['adddipfield'] = v.AdditiveDipoleField
            if hasattr(v, 'DipoleLayerPosition'):
                pars['position'] = v.DipoleLayerPosition
            
        else:
            pars = False
        nc.close()
        return pars
        
    def get_pw(self):
        'return the planewave cutoff used'
        
        ncf = netCDF(self.nc, 'r')
        if 'PlaneWaveCutoff' in ncf.variables:
            pw = ncf.variables['PlaneWaveCutoff'].getValue()
        else:
            pw = None
        ncf.close()
        
        if (isinstance(pw, int)
            or isinstance(pw, float)
            or isinstance(pw,np.int32)):
            return pw
        elif pw is None:
            return None
        else:
            return pw[0]

    def get_dw(self):
        'return the density wave cutoff'
        
        ncf = netCDF(self.nc, 'r')
        if 'Density_WaveCutoff' in ncf.variables:
            dw = ncf.variables['Density_WaveCutoff'].getValue()
        else:
            dw = None
        ncf.close()

        #some old calculations apparently store ints, while newer ones
        #are lists
        if (isinstance(dw, int)
            or isinstance(dw, float)
            or isinstance(dw, np.int32)):
            return dw
        else:
            if dw is None:
                return None
            else:
                return dw[0]
    
    def get_xc(self):
        '''return the self-consistent exchange-correlation functional used

        returns a string'''
        
        nc = netCDF(self.nc, 'r')
        v = 'ExcFunctional'
        if v in nc.variables:
            xc = nc.variables[v][:].tostring().strip()
        else:
            xc = None

        nc.close()
        return xc

    def get_number_of_iterations(self):

        niter = None

        if self.calculation_required():
            self.calculate()

        txt = self.get_txt()
        if os.path.exists(txt):
            f = open(txt,'r')
            for line in f:
                if 'Number of iterations =' in line:
                    niter = int(line.split('=')[1])
                    break
            f.close()
                    
        return niter

    def get_potential_energy(self,
                             atoms=None,
                             force_consistent=False):
        '''
        return the potential energy.
        '''
        
        if self.calculation_required(atoms):
            log.debug('calculation required for energy')
            self.calculate()
        else:
            log.debug('no calculation required for energy')
                        
        nc = netCDF(self.get_nc(), 'r')
        try:
            if force_consistent:
                e = nc.variables['TotalFreeEnergy'][-1]
            else:
                e = nc.variables['TotalEnergy'][-1]
            nc.close()
            return e 
        except (TypeError, KeyError):
            raise RuntimeError('Error in calculating the total energy\n'
                               + 'check %s for error messages'
                               % self.get_txt())

    def get_forces(self, atoms=None):
        """Calculate atomic forces"""
        
        if atoms is None:
            atoms = self.atoms
        if self.calculation_required(atoms):
            self.calculate()
        nc = netCDF(self.get_nc(), 'r')
        forces = nc.variables['DynamicAtomForces'][-1]
        nc.close()
        return forces

    def get_atoms(self):
        'return the atoms attached to a calculator()'
        
        if hasattr(self, 'atoms'):
            if self.atoms is None:
                return None
            atoms = self.atoms.copy()
            #it is not obvious the copy of atoms should have teh same
            #calculator
            atoms.set_calculator(self)
        else:
            atoms = None
        return atoms

    def get_nc(self):
        'return the ncfile used for output'
        
        return self.nc

    def get_txt(self):
        'return the txt file used for output'
        
        if hasattr(self,'txt'):
            return self.txt
        else:
            return None
        
    def get_psp(self, sym=None, z=None):
        '''get the pseudopotential filename from the psp database

        :Parameters:

          sym : string
            the chemical symbol of the species

          z : integer
            the atomic number of the species

        
        you can only specify sym or z. Returns the pseudopotential
        filename, not the full path.
        '''
        if sym is None and z is None:
            return None
        
        if (sym is None and z is not None):
            from ase.data import chemical_symbols
            sym = chemical_symbols[z]
        elif (sym is not None and z is None):
            pass
        else:
            raise Exception, 'You can only specify Z or sym!'
        psp = self.psp[sym]
        return psp
    
    def get_spin_polarized(self):
        'Return True if calculate is spin-polarized or False if not'
        
        #self.calculate() #causes recursion error with get_magnetic_moments
        nc = netCDF(self.nc, 'r')
        if 'ElectronicBands' in nc.variables:
            v = nc.variables['ElectronicBands']
            if hasattr(v, 'SpinPolarization'):
                if v.SpinPolarization == 1:
                    spinpol = False
                elif v.SpinPolarization == 2:
                    spinpol = True
            else:
                spinpol = False
        else:
            spinpol = 'Not defined'

        nc.close()
        return spinpol

    def get_magnetic_moments(self, atoms=None):
        '''return magnetic moments on each atom after the calculation is
        run'''

        if self.calculation_required(atoms):
            self.calculate()
        nc = netCDF(self.nc, 'r')
        if 'InitialAtomicMagneticMoment' in nc.variables:
            mom = nc.variables['InitialAtomicMagneticMoment'][:]
        else:
            mom = [0.0]*len(self.atoms)

        nc.close()
        return mom

    def get_status(self):
        '''get status of calculation from ncfile. usually one of:
        'new',
        'aborted'
        'running'
        'finished'
        None
        '''
        
        nc = netCDF(self.nc, 'r')
        if hasattr(nc, 'status'):
            status = nc.status
        else:
            status = None
        nc.close()
        return status

    def get_calculate_stress(self):
        'return whether stress is calculated or not'
        
        nc = netCDF(self.get_nc(), 'r')
        if 'TotalStress' in nc.variables:
            calcstress = True
        else:
            calcstress = False
        nc.close()
        return calcstress
        
    def get_stress(self, atoms=None):
        '''get stress on the atoms.

        you should have set up the calculation
        to calculate stress first.

        returns [sxx, syy, szz, syz, sxz, sxy]'''
        
        if self.calculation_required(atoms):
            self.calculate()

        nc = netCDF(self.get_nc(), 'r')
        if 'TotalStress' in nc.variables:
            stress = nc.variables['TotalStress'][:]
            #ase expects the 6-element form
            stress = np.take(stress.ravel(), [0, 4, 8, 5, 2, 1])
        else:
            #stress will not be here if you did not set it up by
            #calling set_stress() or in the __init__
            stress = None
        
        nc.close()
        
        return stress

    def get_psp_valence(self, psp):
        '''
        get the psp valence charge on an atom from the pspfile.
        '''
        
        from struct import unpack
        dacapopath = get_dacapopath()

        if os.path.exists(psp):
            #the pspfile may be in the current directory
            #or defined by an absolute path
            fullpsp = psp
        else:
            #or, it is in the default psp path
            fullpsp = os.path.join(dacapopath, psp)

        if os.path.exists(fullpsp.strip()):
            f = open(fullpsp)
            # read past version numbers and text information
            buf = f.read(64)
            # read number valence electrons
            buf = f.read(8)
            fmt = ">d"
            nvalence = unpack(fmt, buf)[0]
            f.close()

        else:
            raise Exception, "%s does not exist" % fullpsp

        return nvalence

    def get_psp_nuclear_charge(self, psp):
        '''
        get the nuclear charge of the atom from teh psp-file.

        This is not the same as the atomic number, nor is it
        necessarily the negative of the number of valence electrons,
        since a psp may be an ion. this function is needed to compute
        centers of ion charge for the dipole moment calculation.

        We read in the valence ion configuration from the psp file and
        add up the charges in each shell.
        '''
        
        from struct import unpack
        dacapopath = get_dacapopath()

        if os.path.exists(psp):
            #the pspfile may be in the current directory
            #or defined by an absolute path
            fullpsp = psp

        else:
            #or, it is in the default psp path
            fullpsp = os.path.join(dacapopath, psp)

        if os.path.exists(fullpsp.strip()):
            f = open(fullpsp)
            unpack('>i', f.read(4))[0]
            for i in range(3):
                f.read(4)
            for i in range(3):
                f.read(4)
            f.read(8)
            f.read(20)
            f.read(8)
            f.read(8)
            f.read(8)
            nvalps = unpack('>i', f.read(4))[0]
            f.read(4)
            f.read(8)
            f.read(8)
            wwnlps = []
            for i in range(nvalps):
                f.read(4)
                wwnlps.append(unpack('>d', f.read(8))[0])
                f.read(8)
            f.close()

        else:
            raise Exception, "%s does not exist" % fullpsp

        return np.array(wwnlps).sum()
       
    def get_valence(self, atoms=None):
        '''return the total number of valence electrons for the
        atoms. valence electrons are read directly from the
        pseudopotentials.

        the psp filenames are stored in the ncfile. They may be just
        the name of the file, in which case the psp may exist in the
        same directory as the ncfile, or in $DACAPOPATH, or the psp
        may be defined by an absolute or relative path. This function
        deals with all these possibilities.
        '''
        
        from struct import unpack
        
        #do not use get_atoms() or recursion occurs
        if atoms is None:
            if hasattr(self, 'atoms'):
                atoms = self.atoms
            else:
                return None

        dacapopath = get_dacapopath()
        totval = 0.0
        for sym in atoms.get_chemical_symbols():
            psp = self.get_psp(sym)
            
            if os.path.exists(psp):
                #the pspfile may be in the current directory
                #or defined by an absolute path
                fullpsp = psp

            #let's also see if we can construct an absolute path to a
            #local or relative path psp.
            abs_path_to_nc = os.path.abspath(self.get_nc())
            base = os.path.split(abs_path_to_nc)[0]
            possible_path_to_psp = os.path.join(base, psp)
            if os.path.exists(possible_path_to_psp):
                fullpsp = possible_path_to_psp
            else:
                #or, it is in the default psp path
                fullpsp = os.path.join(dacapopath, psp)
            if os.path.exists(fullpsp.strip()):
                f = open(fullpsp)
                # read past version numbers and text information
                buf = f.read(64)
                # read number valence electrons
                buf = f.read(8)
                fmt = ">d"
                nvalence = unpack(fmt, buf)[0]
                f.close()
                totval += float(nvalence)
            else:
                print "%s does not exist" % fullpsp
                totval = None
            
        return totval 

    def calculation_required(self, atoms=None, quantities=None):
        '''
        determines if a calculation is needed.

        return True if a calculation is needed to get up to date data.
        return False if no calculation is needed.

        quantities is here because of the ase interface.
        '''
        
        # first, compare if the atoms is the same as the stored atoms
        # if anything has changed, we need to run a calculation
        log.debug('running calculation_required')

        if self.nc is None:
            raise Exception, 'No output ncfile specified!'
                
        if atoms is not None:
            if not self.atoms_are_equal(atoms):
                log.debug('found that atoms != self.atoms')
                tol = 1.0e-6 #tolerance that the unit cell is the same
                new = atoms.get_cell()
                old = self.atoms.get_cell()
                #float comparison of equality
                if not np.all(abs(old-new) < tol): 
                    #this often changes the number of planewaves
                    #which requires a complete restart
                    log.debug('restart required! because cell changed')
                    self.restart()
                else:
                    log.debug('Unitcells apparently the same')
                    
                self.set_atoms(atoms) #we have to update the atoms in any case
                return True
            
        #if we make it past the atoms check, we look in the
        #nc file. if parameters have been changed the status
        #will tell us if a calculation is needed

        #past this point, atoms was None or equal, so there is nothing to
        #update in the calculator

        log.debug('atoms tested equal')
        if os.path.exists(self.nc):
            nc = netCDF(self.nc, 'r')
            if hasattr(nc, 'status'):
                if nc.status == 'finished' and self.ready:
                    nc.close()
                    return False
                elif nc.status == 'running':
                    nc.close()
                    raise DacapoRunning('Dacapo is Running')
                elif nc.status == 'aborted':
                    nc.close()
                    raise DacapoAborted('Dacapo aborted. see txt file!')
                else:
                    log.debug('ncfile exists, but is not ready')
                    nc.close()
                    return True
            else:
                #legacy calculations do not have a status flag in them.
                #let us guess that if the TotalEnergy is there
                #no calculation needs to be run?
                if 'TotalEnergy' in nc.variables:
                    runflag = False
                else:
                    runflag = True
                nc.close()
                log.debug('Legacy calculation')
                return runflag #if no status run calculation
            nc.close()
            
        #default, a calculation is required
        return True

    def get_scratch(self):
        '''finds an appropriate scratch directory for the calculation'''
        
        import getpass
        username = getpass.getuser()

        scratch_dirs = []
        if os.environ.has_key('SCRATCH'):
            scratch_dirs.append(os.environ['SCRATCH'])
        if os.environ.has_key('SCR'):
            scratch_dirs.append(os.environ['SCR'])
        scratch_dirs.append('/scratch/'+username)
        scratch_dirs.append('/scratch/')
        scratch_dirs.append(os.curdir)
        for scratch_dir in scratch_dirs:
            if os.access(scratch_dir, os.W_OK):
                return scratch_dir
        raise IOError, "No suitable scratch directory and no write access \
        to current dir."

    def set_parent(self,parent):
        if hasattr(self,'children'):
            raise RuntimeError,"Cannot create grandparents."
        self.parent = parent

    def attach_child(self,child):
        if hasattr(self,'parent'):
            raise RuntimeError,"Cannot create grandchildren!"
        if not hasattr(self,'children'):
            self.children = []
        self.children.append(child)
        child.set_parent(self)

    def calculate(self):
        '''run a calculation.

        you have to be a little careful with code in here. Use the
        calculation_required function to tell if a calculation is
        required. It is assumed here that if you call this, you mean
        it.'''

        #provide a way to make no calculation get run
        if os.environ.get('DACAPO_DRYRUN', None) is not None:
            raise DacapoDryrun, '$DACAPO_DRYRUN detected, and a calculation \
            attempted'

        if hasattr(self,'children'):
                # We are a parent and call execute_parent_calculation
                self.execute_parent_calculation()
                return

        if hasattr(self,'parent'):  # we're a child and call the parent
                log.debug("I'm a child. Calling parent instead.")
                self.parent.calculate()   # call the parent process to calculate all images
                return
    
        # hack: use the default psp path (see validate.get_dacapopath)
        # export DACAPOPATH to the environment
        env = os.environ
        env['DACAPOPATH'] = get_dacapopath()

        if not self.ready:
            log.debug('Calculator is not ready.')
            if not os.path.exists(self.get_nc()):
                self.initnc()

            log.debug('writing atoms out')
            log.debug(self.atoms)
            self.write_nc() #write atoms to ncfile

            log.debug('writing input out')
            self.write_input() #make sure input is uptodate
   
            #check that the bands get set
            if self.get_nbands() is None:  
                nelectrons = self.get_valence()
                nbands = int(nelectrons * 0.65 + 4)
                self.set_nbands(nbands) 

        log.debug('running a calculation')

        nc = self.get_nc()
        txt = self.get_txt()
        scratch = self.get_scratch()

        if self.stay_alive:
            self.execute_external_dynamics(nc, txt)
            self.ready = True
            self.set_status('finished')
        else:
            # if Dynamics:ExternalIonMotion_script is set in the .nc file from a previous run
            # and stay_alive is false for the continuation run, the Fortran executable continues
            # taking steps of size 0 and ends in an infinite loop.
            # Solution: remove the Dynamics variable if present when not running with stay_alive
            # 
            self.delete_ncattdimvar(self.nc,ncvars=['Dynamics'])
            cmd = "dacapo.run '%(innc)s' -out '%(txt)s' -scratch %(scratch)s"
            cmd = cmd % {'innc':nc,
                         'txt':txt,
                         'scratch':scratch}

            log.debug(cmd)
            # using subprocess instead of commands subprocess is more
            # flexible and works better for stay_alive
            self._dacapo = sp.Popen(cmd,
                                stdout=sp.PIPE,
                                stderr=sp.PIPE,
                                shell=True)
            status = self._dacapo.wait()
            [stdout, stderr] = self._dacapo.communicate()
            output = stdout+stderr
    
            if status is 0: #that means it ended fine!
                self.ready = True
                self.set_status('finished')
            else:
                log.debug('Status was not 0')
                log.debug(output)
                self.ready = False
            # directory cleanup has been moved to self.__del__()
            del self._dacapo

            #Sometimes dacapo dies or is killed abnormally, and in this
            #case an exception should be raised to prevent a geometry
            #optimization from continuing for example. The best way to
            #detect this right now is actually to check the end of the
            #text file to make sure it ends with the right line. The
            #line differs if the job was run in parallel or in serial.
            f = open(txt, 'r')
            lines = f.readlines()
            f.close()

            if 'PAR: msexit halting Master' in lines[-1]:
                pass #standard parallel end
            elif ('TIM' in lines[-2]
                  and 'clexit: exiting the program' in lines[-1]):
                pass #standard serial end
            else:
                # text file does not end as expected, print the last
                # 10 lines and raise exception
                log.debug(string.join(lines[-10:-1], ''))
                s = 'Dacapo output txtfile (%s) did not end normally.\n'
                s += ''.join(lines[-10:-1])
                raise DacapoAbnormalTermination(s % txt)

    def execute_parent_calculation(self):
        '''
        Implementation of an extra level of parallelization, where one jacapo calculator spawns several
        dacapo.run processes. This is used for NEBs parallelized over images.
        '''                

        # hack: use the default psp path (see validate.get_dacapopath)
        # export DACAPOPATH to the environment
        env = os.environ
        env['DACAPOPATH'] = get_dacapopath()

        nchildren = len(self.children)
        log.debug("I'm a parent and start a calculation for ",nchildren," children.")
        self._dacapo = nchildren*[None]
        # export the number of children to the environment
        env = os.environ
        env['JACAPO_NIMAGES'] = str(nchildren)
        
        # start a dacapo.run instance for each child
        for i,child in enumerate(self.children):

            nc = child.get_nc()
            txt= child.get_txt()
            scratch = child.get_scratch()

            if not os.path.exists(nc):
                child.initnc()
            child.write_nc() #write atoms to ncfile
            child.write_input() #make sure input is uptodate

            #check that the bands get set
            if child.get_nbands() is None:
                nelectrons = child.get_valence()
                nbands = int(nelectrons * 0.65 + 4)
                child.set_nbands(nbands)

            env['JACAPO_IMAGE'] = str(i)
            cmd = "dacapo.run  '%(innc)s' -out '%(txt)s' -scratch %(scratch)s"
            cmd = cmd % {'innc':nc,
                         'txt':txt,
                         'scratch':scratch}

            log.debug(cmd)
            self._dacapo[i] = sp.Popen(cmd,stdout=sp.PIPE,stderr=sp.PIPE,shell=True,env=env)
        
        print 'now waiting for all children to finish'
        # now wait for all processes to finish
        for i,child in enumerate(self.children):
            status = self._dacapo[i].wait()
            [stdout,stderr] = self._dacapo[i].communicate()
            output = stdout+stderr
            if status is 0: #that means it ended fine!
                child.ready = True
                child.set_status('finished')
            else:
                log.debug('Status was not 0')
                log.debug(output)
                child.ready = False

        # could also check the end of the output .txt file to make sure everything was fine.
        
        del self._dacapo

    def execute_external_dynamics(self,
                                  nc=None,
                                  txt=None,
                                  stoppfile='stop',
                                  stopprogram=None):
        '''
        Implementation of the stay alive functionality with socket
        communication between dacapo and python.  Known limitations:
        It is not possible to start 2 independent Dacapo calculators
        from the same python process, since the python PID is used as
        identifier for the script[PID].py file.
        '''
        
        from socket import socket, AF_INET, SOCK_STREAM, timeout
        import tempfile
        
        if hasattr(self, "_dacapo"):
            msg = "Starting External Dynamics while Dacapo is runnning: %s"
            msg = msg % str(self._dacapo.poll())
            log.debug(msg)
        else:
            log.debug("No dacapo instance has been started yet")
        log.debug("Stopprogram: %s" % stopprogram)

        if not nc:
            nc = self.get_nc()
        if not txt:
            txt = self.get_txt()
        tempfile.tempdir = os.curdir

        if stopprogram:
            # write stop file
            stfile = open(stoppfile, 'w')
            stfile.write('1 \n')
            stfile.close()

            # signal to dacapo that positions are ready
            # let dacapo continue, it is up to the python mainloop 
            # to allow dacapo enough time to finish properly.
            self._client.send('ok too proceed')

            # Wait for dacapo to acknowledge that netcdf file has
            # been updated, and analysis part of the code has been
            # terminated. Dacapo sends a signal at the end of call
            # clexit().
            log.info("waiting for dacapo to exit...")
            self.s.settimeout(1200.0) # if dacapo exits with an
                                       # error, self.s.accept()
                                       # should time out,
                                      # but we need to give it
                                      # enough time to write the
                                      # wave function to the nc
                                      # file.
            try:
                self._client, self._addr = self.s.accept() # Last
                                                          # mumble
                                                          # before
                                                          # Dacapo
                                                          # dies.
                os.system("sleep 5") # 5 seconds of silence
                                                         # mourning
                                                         # dacapo.
            except timeout:
                print '''Socket connection timed out.'''
                print '''This usually means Dacapo crashed.'''
                
            # close the socket s
            self.s.close()
            self._client.close()

            # remove the script???? file
            ncfile = netCDF(nc, 'r')
            vdyn = ncfile.variables['Dynamics']
            os.system("rm -f '"+vdyn.ExternalIonMotion_script+"'")
            ncfile.close()
            os.system('rm -f '+stoppfile)

            if self._dacapo.poll()==None:  # dacapo is still not dead!
                # but this should do it!
                sp.Popen("kill -9 "+str(self._dacapo.pid), shell=True)
                #if Dacapo dies for example because of too few
                #bands, subprocess never returns an exitcode.
                #very strange, but at least the program is
                #terminated.  print self._dacapo.returncode
            del self._dacapo
            return

        if hasattr(self, '_dacapo') and self._dacapo.poll()==None: 
            # returns  None if  dacapo is running self._dacapo_is_running:

            # calculation_required already updated the positions in
            # the nc file
            self._client.send('ok too proceed')

        else:

            # get process pid that will be used as communication
            # channel
            pid = os.getpid()

            # setup communication channel to dacapo
            from sys    import version
            from string import split
            effpid = (pid)%(2**16-1025)+1025 # This translate pid
                                               # [0;99999] to a number
                                               # in [1025;65535] (the
                                               # allowed socket
                                               # numbers)

            self.s = socket(AF_INET, SOCK_STREAM)
            foundafreesocket = 0
            while not foundafreesocket:
                try:
                    if split(version)[0] > "2":     # new interface
                        self.s.bind(("", effpid))
                    else:                           # old interface
                        self.s.bind("", effpid)
                    foundafreesocket = 1
                except:
                    effpid = effpid + 1

            # write script file that will be used by dacapo
            scriptname = 'script%s.py' % str(pid)
            scriptfile = open(scriptname, 'w')
            scriptfile.write(
"""#!/usr/bin/env python
from socket import *
from sys    import version
from string import split  
s = socket(AF_INET,SOCK_STREAM)
# tell python that dacapo has finished
if split(version)[0] > "2":     # new interface 
     s.connect(("",%(effpid)s))
else:                           # old interface
     s.connect("",%(effpid)s)
# wait for python main loop
s.recv(14)
""" % {'effpid':str(effpid)})
            scriptfile.close()
            os.system('chmod +x ' + scriptname)

            # hack: use the default psp path (see validate.get_dacapopath)
            # export DACAPOPATH to the environment
            env = os.environ
            env['DACAPOPATH'] = get_dacapopath()

            # setup dynamics as external and set the script name
            ncfile = netCDF(nc, 'a')
            if 'Dynamics' not in ncfile.variables:
                vdyn = ncfile.createVariable('Dynamics', 'c', ())
            else:
                vdyn = ncfile.variables['Dynamics']
            vdyn.Type = "ExternalIonMotion" 
            vdyn.ExternalIonMotion_script = './'+ scriptname
            ncfile.close()

            # dacapo is not running start dacapo non blocking
            scratch_in_nc = tempfile.mktemp()
            os.system('mv '+nc+' '+scratch_in_nc)
            os.system('rm -f '+stoppfile)
            scratch = self.get_scratch()
            cmd = "dacapo.run"
            cmd += " '%(innc)s' '%(outnc)s' -out '%(txt)s' -scratch %(scratch)s"
            cmd = cmd % {'innc':scratch_in_nc,
                         'outnc':nc,
                         'txt':txt,
                         'scratch':scratch}

            log.debug(cmd)
            self._dacapo = sp.Popen(cmd,
                                    stdout=sp.PIPE,
                                    stderr=sp.PIPE,
                                    shell=True)

            self.s.listen(1)

        # wait for dacapo  
        self._client, self._addr = self.s.accept()

    def write_nc(self, nc=None, atoms=None):
        '''
        write out atoms to a netcdffile.

        This does not write out the calculation parameters!

        :Parameters:

          nc : string
            ncfilename to write to. this file will get clobbered
            if it already exists.

          atoms : ASE.Atoms
            atoms to write. if None use the attached atoms
            if no atoms are attached only the calculator is
            written out. 

        the ncfile is always opened in 'a' mode.

        note: it is good practice to use the atoms argument to make
        sure that the geometry you mean gets written! Otherwise, the
        atoms in the calculator is used, which may be different than
        the external copy of the atoms.

        '''

        log.debug('writing atoms to ncfile with write_nc')
        #no filename was provided to function, use the current ncfile
        if nc is None: 
            nc = self.get_nc()
        
        if nc != self.nc:
            #this means we are writing a new file, and we should copy
            #the old file to it first.  this makes sure the old
            #calculator settings are preserved
            new = nc
            old = self.nc
            log.debug('Copying old ncfile to new ncfile')
            log.debug("cp '%s' '%s'" % (old, new))
            os.system("cp '%s' '%s'" % (old, new))
              
        if atoms is None:
            atoms = self.get_atoms()

        log.debug('self.atoms = %s' % str(self.atoms))
        log.debug('atoms = %s' % str(atoms))

        if atoms is not None: #there may still be no atoms attached
            log.debug('about to write to %s' % nc)
            ncf = netCDF(nc, 'a')

            if 'number_of_dynamic_atoms' not in ncf.dimensions:
                ncf.createDimension('number_of_dynamic_atoms',
                                    len(atoms))
            else:
                # number of atoms is already a dimension, but we might
                # be setting new atoms here
                # check for same atom symbols (implicitly includes
                # a length check)
                symbols = np.array(['%2s' % s for s in
                                    atoms.get_chemical_symbols()], dtype='c')
                ncsym = ncf.variables['DynamicAtomSpecies'][:]
                if (symbols.size != ncsym.size) or (np.any(ncsym != symbols)):
                    # the number of atoms or their order has changed.
                    # Treat this as a new calculation and reset
                    # number_of_ionic_steps and
                    # number_of_dynamic_atoms.
                    ncf.close() #nc file must be closed for
                                 #delete_ncattdimvar to work correctly
                    self.delete_ncattdimvar(nc, ncattrs=[],
                                            ncdims=['number_of_dynamic_atoms',
                                                    'number_ionic_steps'])
                    ncf = netCDF(nc, 'a')
                    ncf.createDimension('number_of_dynamic_atoms',
                                                  len(atoms)) 
                    ncf.createDimension('number_ionic_steps', None)
                    self._set_frame_number(0)
                    ncf.close() #nc file must be closed for restart to
                                #work correctly
                    self.restart()
                    ncf = netCDF(nc, 'a')

            #now, create variables
            if 'DynamicAtomSpecies' not in ncf.variables:
                sym = ncf.createVariable('DynamicAtomSpecies',
                                         'c',
                                         ('number_of_dynamic_atoms',
                                          'dim2',))
            else:
                sym = ncf.variables['DynamicAtomSpecies']

            #note explicit array casting was required here
            symbols = atoms.get_chemical_symbols()
            sym[:] = np.array(['%2s' % s for s in symbols], dtype='c')

            if 'DynamicAtomPositions' not in ncf.variables:
                pos = ncf.createVariable('DynamicAtomPositions',
                                         'd',
                                         ('number_ionic_steps',
                                          'number_of_dynamic_atoms',
                                          'dim3'))
            else:
                pos = ncf.variables['DynamicAtomPositions']

            spos = atoms.get_scaled_positions()
            if pos.typecode() == 'f':
                spos = np.array(spos, dtype=np.float32)
            pos[self._frame, :] = spos

            if 'UnitCell' not in ncf.variables:
                uc = ncf.createVariable('UnitCell', 'd',
                                        ('number_ionic_steps',
                                         'dim3', 'dim3'))
            else:
                uc = ncf.variables['UnitCell']

            cell = atoms.get_cell()
            if uc.typecode() == 'f':
                cell = np.array(cell, dtype=np.float32)
                
            uc[self._frame, :] = cell

            if 'AtomTags' not in ncf.variables:
                tags = ncf.createVariable('AtomTags', 'i',
                                          ('number_of_dynamic_atoms',))
            else:
                tags = ncf.variables['AtomTags']

            tags[:] = np.array(atoms.get_tags(), np.int32)

            if 'InitialAtomicMagneticMoment' not in ncf.variables:
                mom = ncf.createVariable('InitialAtomicMagneticMoment',
                                         'd',
                                         ('number_of_dynamic_atoms',))
            else:
                mom = ncf.variables['InitialAtomicMagneticMoment']

            #explain why we have to use get_initial_magnetic_moments()
            moms = atoms.get_initial_magnetic_moments()
            if mom.typecode() == 'f':
                moms = np.array(moms, dtype=np.float32)
            mom[:] = moms
            
            #finally the atom pseudopotentials
            for sym in atoms.get_chemical_symbols():
                vn = 'AtomProperty_%s' % sym
                if vn not in ncf.variables:
                    p = ncf.createVariable(vn, 'c', ('dim20',))
                else:
                    p = ncf.variables[vn]

                ppath = self.get_psp(sym=sym)
                p.PspotFile = ppath

            ncf.sync()
            ncf.close()
               
            #store constraints if they exist
            constraints = atoms._get_constraints()
            if constraints != []:
                nc = netCDF(self.get_nc(), 'a')
                if 'constraints' not in nc.variables:
                    if 'dim1' not in nc.dimensions:
                        nc.createDimension('dim1', 1)
                    c = nc.createVariable('constraints', 'c', ('dim1',))
                else:
                    c = nc.variables['constraints']
                #we store the pickle string as an attribute of a
                #netcdf variable because that way we do not have to
                #know how long the string is.  with a character
                #variable you have to specify the dimension of the
                #string ahead of time.
                c.data = pickle.dumps(constraints)
                nc.close()
            else:
                # getting here means there where no constraints on the
                # atoms just written we should check if there are any
                # old constraints left in the ncfile
                # from a previous atoms, and delete them if so
                delete_constraints = False
                nc = netCDF(self.get_nc())
                if 'constraints' in nc.variables:
                    delete_constraints = True
                nc.close()

                if delete_constraints:
                    log.debug('deleting old constraints')
                    self.delete_ncattdimvar(self.nc,
                                            ncvars=['constraints'])
        
    def read_atoms(filename):
        '''read atoms and calculator from an existing netcdf file.

        :Parameters:

          filename : string
            name of file to read from.

        static method

        example::
        
          >>> atoms = Jacapo.read_atoms(ncfile)
          >>> calc = atoms.get_calculator()

        this method is here for legacy purposes. I used to use it alot.
        '''
        
        calc = Jacapo(filename)
        atoms = calc.get_atoms()
        return atoms
    
    read_atoms = staticmethod(read_atoms)

    def read_only_atoms(self, ncfile):
        '''read only the atoms from an existing netcdf file. Used to
        initialize a calculator from a ncfilename.

        :Parameters:

          ncfile : string
            name of file to read from.

        return ASE.Atoms with no calculator attached or None if no
        atoms found
        '''
        
        from ase import Atoms
        
        nc = netCDF(ncfile, 'r')
        #some ncfiles do not have atoms in them
        if 'UnitCell' not in nc.variables:
            log.debug('no unit cell found in ncfile')
            nc.close()
            return None
        
        cell = nc.variables['UnitCell'][:][-1]
        sym = nc.variables['DynamicAtomSpecies'][:]
        symbols = [x.tostring().strip() for x in sym]
        spos = nc.variables['DynamicAtomPositions'][:][-1]

        pos = np.dot(spos, cell)
        
        atoms = Atoms(symbols=symbols,
                      positions=pos,
                      cell=cell)

        if 'AtomTags' in nc.variables:
            tags = nc.variables['AtomTags'][:]
            atoms.set_tags(tags)

        if 'InitialAtomicMagneticMoment' in nc.variables:
            mom = nc.variables['InitialAtomicMagneticMoment'][:]
            atoms.set_initial_magnetic_moments(mom)

        #update psp database
        for sym in symbols:
            vn = 'AtomProperty_%s' % sym
            if vn in nc.variables:
                var = nc.variables[vn]
                pspfile = var.PspotFile
                self.psp[sym] = pspfile

        #get constraints if they exist
        c = nc.variables.get('constraints', None)
        if c is not None:
            constraints = pickle.loads(c.data)
            atoms.set_constraint(constraints)
                    
        nc.close()
        
        return atoms
        
    def delete_ncattdimvar(self, ncf, ncattrs=None, ncdims=None, ncvars=None):
        '''
        helper function to delete attributes,
        dimensions and variables in a netcdffile

        this functionality is not implemented for some reason in
        netcdf, so the only way to do this is to copy all the
        attributes, dimensions, and variables to a new file, excluding
        the ones you want to delete and then rename the new file.

        if you delete a dimension, all variables with that dimension
        are also deleted.
        '''

        if ncattrs is None:
            ncattrs = []
        if ncdims is None:
            ncdims = []
        if ncvars is None:
            ncvars = []
        
        log.debug('beginning: going to delete dims: %s' % str(ncdims))
        log.debug('beginning: going to delete vars: %s' % str(ncvars))
            
        oldnc = netCDF(ncf, 'r')

        #h,tempnc = tempfile.mkstemp(dir='.',suffix='.nc')
        tempnc = ncf+'.temp'
        
        newnc = netCDF(tempnc, 'w')

        for attr in dir(oldnc):
            if attr in ['close', 'createDimension',
                        'createVariable', 'flush', 'sync']:
                continue
            if attr in ncattrs:
                continue #do not copy this attribute
            setattr(newnc, attr, getattr(oldnc, attr))
           
        #copy dimensions
        for dim in oldnc.dimensions:
            if dim in ncdims:
                log.debug('deleting %s of %s' % (dim, str(ncdims)))
                continue #do not copy this dimension
            size = oldnc.dimensions[dim]
            
            newnc.createDimension(dim, size)

        # we need to delete all variables that depended on a deleted dimension
        for v in oldnc.variables:
            dims1 = oldnc.variables[v].dimensions
            for dim in ncdims:
                if dim in dims1:
                    s = 'deleting "%s" because it depends on dim "%s"'
                    log.debug(s %(v, dim))
                    ncvars.append(v)

        #copy variables, except the ones to delete
        for v in oldnc.variables:
            if v in ncvars:
                log.debug('vars to delete: %s ' % ncvars)
                log.debug('deleting ncvar: %s' % v)
                continue #we do not copy this v over

            ncvar = oldnc.variables[v]
            tcode = ncvar.typecode()
            #char typecodes do not come out right apparently
            if tcode == " ":
                tcode = 'c'
            
            ncvar2 = newnc.createVariable(v, tcode, ncvar.dimensions)
            try:
                ncvar2[:] = ncvar[:]
            except TypeError:
                #this exception occurs for scalar variables
                #use getValue and assignValue instead
                ncvar2.assignValue(ncvar.getValue())

            #and variable attributes
            #print dir(ncvar)
            for att in dir(ncvar):
                if att in ['assignValue', 'getValue', 'typecode']:
                    continue
                setattr(ncvar2, att, getattr(ncvar, att))

        oldnc.close()
        newnc.close()

        s = 'looking for .nfs files before copying: %s'
        log.debug(s % glob.glob('.nfs*'))
        
        #ack!!! this makes .nfsxxx files!!!
        #os.close(h) #this avoids the stupid .nfsxxx file
        #import shutil
        #shutil.move(tempnc,ncf)

        #this seems to avoid making the .nfs files 
        os.system("cp '%s' '%s'" % (tempnc, ncf))
        os.system("rm '%s'" % tempnc)

        s = 'looking for .nfs files after copying: %s'
        log.debug(s %  glob.glob('.nfs*'))
        
    def restart(self):
        '''
        Restart the calculator by deleting nc dimensions that will
        be rewritten on the next calculation. This is sometimes required
        when certain dimensions change related to unitcell size changes
        planewave/densitywave cutoffs and kpt changes. These can cause
        fortran netcdf errors if the data does not match the pre-defined
        dimension sizes.

        also delete all the output from previous calculation.
        '''
        
        log.debug('restarting!')
        
        ncdims = ['number_plane_waves',
                  'number_IBZ_kpoints',
                  'softgrid_dim1',
                  'softgrid_dim2',
                  'softgrid_dim3',
                  'hardgrid_dim1',
                  'hardgrid_dim2',
                  'hardgrid_dim3',
                  'max_projectors_per_atom',
                  'atomdos_energygrid_size',
                  'atomdos_angular_channels',
                  'atomdos_radial_orbs']

        ncvars = ['TotalEnergy',
                  'TotalFreeEnergy',
                  'EvaluateTotalEnergy',
                  'DynamicAtomForces',
                  'FermiLevel',
                  'EnsembleXCEnergies',
                  'AtomProjectedDOS_IntegratedDOS',
                  'AtomProjectedDOS_OrdinalMap',
                  'NumberPlaneWavesKpoint',
                  'AtomProjectedDOS_EnergyResolvedDOS',
                  'AtomProjectedDOS_EnergyGrid',
                  'EvaluateCorrelationEnergy',
                  'DynamicAtomVelocities',
                  'KpointWeight',
                  'EvaluateExchangeEnergy',
                  'EffectivePotential',
                  'TotalStress',
                  'ChargeDensity',
                  'WaveFunction',
                  'WaveFunctionFFTindex',
                  'NumberOfNLProjectors',
                  'NLProjectorPsi',
                  'TypeNLProjector1',
                  'NumberofNLProjectors',
                  'PartialCoreDensity',
                  'ChargeDensity',
                  'ElectrostaticPotential',
                  'StructureFactor',
                  'EigenValues',
                  'OccupationNumbers']
                  
        self.delete_ncattdimvar(self.nc,
                                ncattrs=[],
                                ncdims=ncdims,
                                ncvars=ncvars)

        self.set_status('new')
        self.ready = False

    def get_convergence(self):
        'return convergence settings for Dacapo'
        
        nc = netCDF(self.get_nc(), 'r')
        vname = 'ConvergenceControl'
        if vname in nc.variables:
            v = nc.variables[vname]
            convergence = {}
            if hasattr(v, 'AbsoluteEnergyConvergence'):
                convergence['energy'] = v.AbsoluteEnergyConvergence[0]
            if hasattr(v, 'DensityConvergence'):
                convergence['density'] = v.DensityConvergence[0]
            if hasattr(v, 'OccupationConvergence'):
                convergence['occupation'] = v.OccupationConvergence[0]
            if hasattr(v, 'MaxNumberOfSteps'):
                convergence['maxsteps'] = v.MaxNumberOfSteps[0]
            if hasattr(v, 'CPUTimeLimit'):
                convergence['cputime'] = v.CPUTimeLimit[0]
        else:
            convergence = None

        nc.close()
        return convergence

    def set_convergence(self,
                        energy=0.00001,
                        density=0.0001,
                        occupation=0.001,
                        maxsteps=None,
                        maxtime=None
                        ):
        '''set convergence criteria for stopping the dacapo calculator.

        :Parameters:

          energy : float
            set total energy change (eV) required for stopping

          density : float
            set density change required for stopping

          occupation : float
            set occupation change required for stopping

          maxsteps : integer
            specify maximum number of steps to take

          maxtime : integer
            specify maximum number of hours to run.

        Autopilot not supported here.
        '''
        
        nc = netCDF(self.get_nc(), 'a')
        vname = 'ConvergenceControl'
        if vname in nc.variables:
            v = nc.variables[vname]
        else:
            v = nc.createVariable(vname, 'c', ('dim1',))

        if energy is not None:
            v.AbsoluteEnergyConvergence = energy
        if density is not None:
            v.DensityConvergence = density
        if occupation is not None:
            v.OccupationConvergence = occupation
        if maxsteps is not None:
            v.MaxNumberOfSteps = maxsteps
        if maxtime is not None:
            v.CPUTimeLimit = maxtime

        nc.sync()
        nc.close()

    def get_charge_mixing(self):
        'return charge mixing parameters'
        
        nc = netCDF(self.get_nc(), 'r')
        vname = 'ChargeMixing'
        if vname in nc.variables:
            v = nc.variables[vname]
            charge_mixing = {}
            if hasattr(v, 'Method'):
                charge_mixing['method'] = v.Method
            if hasattr(v, 'UpdateCharge'):
                charge_mixing['updatecharge'] = v.UpdateCharge
            if hasattr(v, 'Pulay_MixingHistory'):
                charge_mixing['mixinghistory'] = v.Pulay_MixingHistory[0]
            if hasattr(v, 'Pulay_DensityMixingCoeff'):
                charge_mixing['mixingcoeff'] = v.Pulay_DensityMixingCoeff[0]
            if hasattr(v, 'Pulay_KerkerPrecondition'):
                charge_mixing['precondition'] = v.Pulay_KerkerPrecondition
        else:
            charge_mixing = None

        nc.close()
        return charge_mixing

    def set_charge_mixing(self,
                          method='Pulay',
                          mixinghistory=10,
                          mixingcoeff=0.1,
                          precondition='No',
                          updatecharge='Yes'):
        '''set density mixing method and parameters

        :Parameters:

          method : string
            'Pulay' for Pulay mixing. only one supported now

          mixinghistory : integer
            number of iterations to mix
            Number of charge residual vectors stored for generating
            the Pulay estimate on the self-consistent charge density,
            see Sec. 4.2 in Kresse/Furthmuller:
            Comp. Mat. Sci. 6 (1996) p34ff

          mixingcoeff : float
            Mixing coefficient for Pulay charge mixing, corresponding
            to A in G$^1$ in Sec. 4.2 in Kresse/Furthmuller:
            Comp. Mat. Sci. 6 (1996) p34ff
                        
          precondition : string
            'Yes' or 'No'
            
            * "Yes" : Kerker preconditiong is used,
               i.e. q$_0$ is different from zero, see eq. 82
               in Kresse/Furthmuller: Comp. Mat. Sci. 6 (1996).
               The value of q$_0$ is fix to give a damping of 20
               of the lowest q vector.
            
            * "No" : q$_0$ is zero and mixing is linear (default).

          updatecharge : string
            'Yes' or 'No'
            
            * "Yes" : Perform charge mixing according to
               ChargeMixing:Method setting
              
            * "No" : Freeze charge to initial value.
               This setting is useful when evaluating the Harris-Foulkes
               density functional
              
        '''
        
        if method == 'Pulay':
            nc = netCDF(self.get_nc(), 'a')
            vname = 'ChargeMixing'
            if vname in nc.variables:
                v = nc.variables[vname]
            else:
                v = nc.createVariable(vname, 'c', ('dim1',))

            v.Method = 'Pulay'
            v.UpdateCharge = updatecharge
            v.Pulay_MixingHistory = mixinghistory
            v.Pulay_DensityMixingCoeff = mixingcoeff
            v.Pulay_KerkerPrecondition = precondition

            nc.sync()
            nc.close()

        self.ready = False

    def set_electronic_minimization(self,
                                    method='eigsolve',
                                    diagsperband=2):
        '''set the eigensolver method

        Selector for which subroutine to use for electronic
        minimization

        Recognized options : "resmin", "eigsolve" and "rmm-diis".

        * "resmin" : Power method (Lennart Bengtson), can only handle
           k-point parallization.

        * "eigsolve : Block Davidson algorithm
           (Claus Bendtsen et al).

        * "rmm-diis : Residual minimization
           method (RMM), using DIIS (direct inversion in the iterate
           subspace) The implementaion follows closely the algorithm
           outlined in Kresse and Furthmuller, Comp. Mat. Sci, III.G/III.H
        
        :Parameters:

          method : string
            should be 'resmin', 'eigsolve' or 'rmm-diis'

          diagsperband : int
            The number of diagonalizations per band for
            electronic minimization algorithms (maps onto internal
            variable ndiapb). Applies for both
            ElectronicMinimization:Method = "resmin" and "eigsolve".
            default value = 2
        '''
        
        nc = netCDF(self.get_nc(), 'a')
        
        vname = 'ElectronicMinimization'
        if vname in nc.variables:
            v = nc.variables[vname]
        else:
            log.debug('Creating ElectronicMinimization')
            v = nc.createVariable(vname, 'c', ('dim1',))

        log.debug('setting method for ElectronicMinimization: % s' % method)
        v.Method = method
        log.debug('setting DiagonalizationsBand for ElectronicMinimization')
        if diagsperband is not None:
            v.DiagonalizationsPerBand = diagsperband

        log.debug('synchronizing ncfile')
        nc.sync()
        
        nc.close()

    def get_electronic_minimization(self):
        '''get method and diagonalizations per band for electronic
        minimization algorithms'''

        log.debug('getting electronic minimization parameters')
        
        nc = netCDF(self.get_nc(), 'r')
        vname = 'ElectronicMinimization'
        if vname in nc.variables:
            v = nc.variables[vname]
            method = v.Method
            if hasattr(v, 'DiagonalizationsPerBand'):
                diagsperband = v.DiagonalizationsPerBand[0]
            else:
                diagsperband = None
        else:
            method = None
            diagsperband = None
        nc.close()
        return {'method':method,
                'diagsperband':diagsperband}

    def get_occupationstatistics(self):
        'return occupation statistics method'
        
        nc = netCDF(self.get_nc(), 'r')
        if 'ElectronicBands' in nc.variables:
            v = nc.variables['ElectronicBands']
            if hasattr(v, 'OccupationStatistics'):
                occstat = v.OccupationStatistics
            else:
                occstat = None
        else:
            occstat = None
        nc.close()
        return occstat

    def set_occupationstatistics(self, method):
        '''
        set the method used for smearing the occupations.

        :Parameters:

          method : string
            one of 'FermiDirac' or 'MethfesselPaxton'
            Currently, the Methfessel-Paxton scheme (PRB 40, 3616 (1989).)
            is implemented to 1th order (which is recommemded by most authors).
            'FermiDirac' is the default
        '''
            
        nc = netCDF(self.get_nc(), 'a')
        if 'ElectronicBands' in nc.variables:
            v = nc.variables['ElectronicBands']
            v.OccupationStatistics = method
        
        nc.sync()
        nc.close()

    def get_fermi_level(self):
        'return Fermi level'
        
        if self.calculation_required():
            self.calculate()
        nc = netCDF(self.get_nc(), 'r')
        ef = nc.variables['FermiLevel'][-1]
        nc.close()
        return ef

    def get_occupation_numbers(self, kpt=0, spin=0):
        '''return occupancies of eigenstates for a kpt and spin

        :Parameters:

          kpt : integer
            index of the IBZ kpoint you want the occupation of

          spin : integer
            0 or 1
        '''
        
        if self.calculation_required():
            self.calculate()
        nc = netCDF(self.get_nc(), 'r')
        occ = nc.variables['OccupationNumbers'][:][-1][kpt, spin]
        nc.close()
        return occ

    def get_xc_energies(self, *functional):
        """
        Get energies for different functionals self-consistent and
        non-self-consistent.

        :Parameters:

          functional : strings
            some set of 'PZ','VWN','PW91','PBE','revPBE', 'RPBE'

        This function returns the self-consistent energy and/or
        energies associated with various functionals. 
        The functionals are currently PZ,VWN,PW91,PBE,revPBE, RPBE.
        The different energies may be useful for calculating improved
        adsorption energies as in B. Hammer, L.B. Hansen and
        J.K. Norskov, Phys. Rev. B 59,7413. 
        Examples: 
        get_xcenergies() #returns all the energies
        get_xcenergies('PBE') # returns the PBE total energy
        get_xcenergies('PW91','PBE','revPBE') # returns a
        # list of energies in the order asked for
        """
        
        if self.calculation_required():
            self.calculate()

        nc = netCDF(self.get_nc(), 'r')

        funcenergies = nc.variables['EvaluateTotalEnergy'][:][-1]
        xcfuncs = nc.variables['EvalFunctionalOfDensity_XC'][:]

        nc.close()
        
        xcfuncs = [xc.tostring().strip() for xc in xcfuncs]
        edict = dict(zip(xcfuncs, funcenergies))

        if len(functional) == 0:
            #get all energies by default
            functional = xcfuncs

        return [edict[xc] for xc in functional]

    # break of compatibility
    def get_ados_data(self,
                      atoms,
                      orbitals,
                      cutoff,
                      spin):
        '''get atom projected data

        :Parameters:

          atoms  
              list of atom indices (integers)

          orbitals
              list of strings
              ['s','p','d'],
              ['px','py','pz']
              ['d_zz', 'dxx-yy', 'd_xy', 'd_xz', 'd_yz']

          cutoff : string
            cutoff radius you want the results for 'short' or 'infinite'

          spin
            : list of integers
            spin you want the results for
            [0] or [1] or [0,1] for both

        returns (egrid, ados)
        egrid has the fermi level at 0 eV
        '''
              
        if self.calculation_required():
            self.calculate()
        nc = netCDF(self.get_nc(), 'r')
        omapvar = nc.variables['AtomProjectedDOS_OrdinalMap']
        omap = omapvar[:] #indices
        c = omapvar.AngularChannels
        channels = [x.strip() for x in c.split(',')] #channel names
        #this has dimensions(nprojections, nspins, npoints)
        ados = nc.variables['AtomProjectedDOS_EnergyResolvedDOS'][:]
        #this is the energy grid for all the atoms
        egrid = nc.variables['AtomProjectedDOS_EnergyGrid'][:]
        nc.close()

        #it is apparently not necessary to normalize the egrid to
        #the Fermi level. the data is already for ef = 0.

        #get list of orbitals, replace 'p' and 'd' in needed
        orbs = []
        for o in orbitals:
            if o == 'p':
                orbs += ['p_x', 'p_y', 'p_z']
            elif o == 'd':
                orbs += ['d_zz', 'dxx-yy', 'd_xy', 'd_xz', 'd_yz']
            else:
                orbs += [o]

        orbinds = [channels.index(x) for x in orbs]

        cutdict = {'infinite':0,
                   'short':1}

        icut = cutdict[cutoff]

        ydata = np.zeros(len(egrid), np.float)
        
        for atomind in atoms:
            for oi in orbinds:
                ind = omap[atomind, icut, oi]

                for si in spin:
                    ydata += ados[ind, si]

        return (egrid, ydata)

    def get_all_eigenvalues(self, spin=0):
        '''return all the eigenvalues at all the kpoints for a spin.

        :Parameters:

          spin : integer
            which spin the eigenvalues are for'''
        
        if self.calculation_required():
            self.calculate()
        nc = netCDF(self.get_nc(), 'r')
        ev = nc.variables['EigenValues'][:][-1][:, spin]
        nc.close()
        return ev
    
    def get_eigenvalues(self, kpt=0, spin=0):
        '''return the eigenvalues for a kpt and spin

        :Parameters:

          kpt : integer
            index of the IBZ kpoint

          spin : integer
            which spin the eigenvalues are for'''
        
        if self.calculation_required():
            self.calculate()
        nc = netCDF(self.get_nc(), 'r')
        ev = nc.variables['EigenValues'][:][-1][kpt, spin]
        nc.close()
        return ev
    
    def get_k_point_weights(self):
        'return the weights on the IBZ kpoints'
        
        if self.calculation_required():
            self.calculate()
        nc = netCDF(self.get_nc(), 'r')
        kw = nc.variables['KpointWeight'][:]
        nc.close()
        return kw

    def get_magnetic_moment(self, atoms=None):
        'calculates the magnetic moment (Bohr-magnetons) of the supercell'

        if not self.get_spin_polarized():
            return None
        
        if self.calculation_required():
            self.calculate()

        nibzk = len(self.get_ibz_kpoints())
        ibzkw = self.get_k_point_weights()
        spinup, spindn = 0.0, 0.0

        for k in range(nibzk):

            spinup += self.get_occupation_numbers(k, 0).sum()*ibzkw[k]
            spindn += self.get_occupation_numbers(k, 1).sum()*ibzkw[k]

        return (spinup - spindn)

    def get_number_of_spins(self):
        'if spin-polarized returns 2, if not returns 1'
        
        if self.calculation_required():
            self.calculate()
        nc = netCDF(self.get_nc(), 'r')
        spv = nc.variables['ElectronicBands']
        nc.close()

        if hasattr(spv, 'SpinPolarization'):
            return spv.SpinPolarization
        else:
            return 1

    def get_ibz_kpoints(self):
        'return list of kpoints in the irreducible brillouin zone'
        
        if self.calculation_required():
            self.calculate()
        nc = netCDF(self.get_nc(), 'r')
        ibz = nc.variables['IBZKpoints'][:]
        nc.close()
        return ibz

    get_ibz_k_points = get_ibz_kpoints

    def get_bz_k_points(self):
        'return list of kpoints in the Brillouin zone'
        
        nc = netCDF(self.get_nc(), 'r')
        if 'BZKpoints' in nc.variables:
            bz = nc.variables['BZKpoints'][:]
        else:
            bz = None
        nc.close()
        return bz
    
    def get_effective_potential(self, spin=1):
        '''
        returns the realspace local effective potential for the spin.
        the units of the potential are eV

        :Parameters:

          spin : integer
             specify which spin you want, 0 or 1
          
        '''
        
        if self.calculation_required():
            self.calculate()
            
        nc = netCDF(self.get_nc(), 'r')
        efp = np.transpose(nc.variables['EffectivePotential'][:][spin])
        nc.close()
        fftgrids = self.get_fftgrid()
        hardgrid = fftgrids['hard']
        x, y, z = self.get_ucgrid(hardgrid)
        return (x, y, z, efp)
        
    def get_electrostatic_potential(self, spin=0):
        '''get electrostatic potential

        Netcdf documentation::
        
          double ElectrostaticPotential(number_of_spin,
                                        hardgrid_dim3,
                                        hardgrid_dim2,
                                        hardgrid_dim1) ;
                 ElectrostaticPotential:
                     Description = "realspace local effective potential" ;
                     unit = "eV" ;
                     
        '''
        
        if self.calculation_required():
            self.calculate()
            
        nc = netCDF(self.get_nc(), 'r')
        esp = np.transpose(nc.variables['ElectrostaticPotential'][:][spin])
        nc.close()
        fftgrids = self.get_fftgrid()

        x, y, z = self.get_ucgrid(fftgrids['hard'])
        
        return (x, y, z, esp)
    
    def get_charge_density(self, spin=0):
        '''
        return x,y,z,charge density data
        
        x,y,z are grids sampling the unit cell
        cd is the charge density data

        netcdf documentation::
        
          ChargeDensity(number_of_spin,
                        hardgrid_dim3,
                        hardgrid_dim2,
                        hardgrid_dim1)
          ChargeDensity:Description = "realspace charge density" ;
                  ChargeDensity:unit = "-e/A^3" ;

        '''
        
        if self.calculation_required():
            self.calculate()
            
        nc = netCDF(self.get_nc(), 'r')
     
        cd = np.transpose(nc.variables['ChargeDensity'][:][spin])

        #I am not completely sure why this has to be done
        #it does give units of electrons/ang**3
        vol = self.get_atoms().get_volume()
        cd /= vol
        nc.close()
        grids = self.get_fftgrid()

        x, y, z = self.get_ucgrid(grids['hard'])
        return x, y, z, cd

    def get_ucgrid(self, dims):
        '''Return X,Y,Z grids for uniform sampling of the unit cell

        dims = (n0,n1,n2)

        n0 points along unitcell vector 0
        n1 points along unitcell vector 1
        n2 points along unitcell vector 2
        '''
        
        n0, n1, n2 = dims
        
        s0 = 1.0/n0
        s1 = 1.0/n1
        s2 = 1.0/n2

        X, Y, Z = np.mgrid[0.0:1.0:s0,
                          0.0:1.0:s1,
                          0.0:1.0:s2]
        
        C = np.column_stack([X.ravel(),
                             Y.ravel(),
                             Z.ravel()])
        
        atoms = self.get_atoms()
        uc = atoms.get_cell()
        real = np.dot(C, uc)

        #now convert arrays back to unitcell shape
        RX = np.reshape(real[:, 0], (n0, n1, n2))
        RY = np.reshape(real[:, 1], (n0, n1, n2))
        RZ = np.reshape(real[:, 2], (n0, n1, n2))
        return (RX, RY, RZ)

    def get_number_of_grid_points(self):
        'return soft fft grid'
        
        # needed by ase.dft.wannier
        fftgrids = self.get_fftgrid()
        return np.array(fftgrids['soft'])
    
    def get_wannier_localization_matrix(self, nbands, dirG, kpoint,
                                        nextkpoint, G_I, spin):
        'return wannier localization  matrix'

        if self.calculation_required():
            self.calculate()

        if not hasattr(self, 'wannier'):
            from utils.wannier import Wannier
            self.wannier = Wannier(self)
            self.wannier.set_bands(nbands)
            self.wannier.set_spin(spin)
        locmat = self.wannier.get_zi_bloch_matrix(dirG,
                                                  kpoint,
                                                  nextkpoint,
                                                  G_I)
        return locmat

    def initial_wannier(self,
                        initialwannier,
                        kpointgrid,
                        fixedstates,
                        edf,
                        spin):
        'return initial wannier'

        if self.calculation_required():
            self.calculate()

        if not hasattr(self, 'wannier'):
            from utils.wannier import Wannier
            self.wannier = Wannier(self)

        self.wannier.set_data(initialwannier)
        self.wannier.set_k_point_grid(kpointgrid)
        self.wannier.set_spin(spin)

        waves = [[self.get_reciprocal_bloch_function(band=band,
                                                     kpt=kpt,
                                                     spin=spin)
                  for band in range(self.get_nbands())]
                  for kpt in range(len(self.get_ibz_k_points()))]

        self.wannier.setup_m_matrix(waves, self.get_bz_k_points())

        #lfn is too keep line length below 78 characters
        lfn = self.wannier.get_list_of_coefficients_and_rotation_matrices
        c, U = lfn((self.get_nbands(), fixedstates, edf))

        U = np.array(U)
        for k in range(len(c)):
            c[k] = np.array(c[k])
        return c, U

    def get_dipole_moment(self,atoms=None):
        '''
        return dipole moment of unit cell

        Defined by the vector connecting the center of electron charge
        density to the center of nuclear charge density.

        Units = eV*angstrom

        1 Debye = 0.208194 eV*angstrom

        '''
        if self.calculation_required():
            self.calculate()

        if atoms is None:
            atoms = self.get_atoms()
        
        #center of electron charge density
        x, y, z, cd = self.get_charge_density()

        n1, n2, n3 = cd.shape
        nelements = n1*n2*n3
        voxel_volume = atoms.get_volume()/nelements
        total_electron_charge = -cd.sum()*voxel_volume

        
        electron_density_center = np.array([(cd*x).sum(),
                                            (cd*y).sum(),
                                            (cd*z).sum()])
        electron_density_center *= voxel_volume
        electron_density_center /= total_electron_charge
       
        electron_dipole_moment = electron_density_center*total_electron_charge
        electron_dipole_moment *= -1.0 #we need the - here so the two
                                        #negatives don't cancel
        # now the ion charge center
        psps = self.get_pseudopotentials()['pspdict']
        ion_charge_center = np.array([0.0, 0.0, 0.0])
        total_ion_charge = 0.0
        for atom in atoms:
            Z = self.get_psp_nuclear_charge(psps[atom.symbol])
            total_ion_charge += Z
            pos = atom.position
            ion_charge_center += Z*pos

        ion_charge_center /= total_ion_charge
        ion_dipole_moment = ion_charge_center*total_ion_charge

        dipole_vector = (ion_dipole_moment + electron_dipole_moment)
        return dipole_vector

    
    def get_reciprocal_bloch_function(self, band=0, kpt=0, spin=0):
        '''return the reciprocal bloch function. Need for Jacapo
        Wannier class.'''
        
        if self.calculation_required():
            self.calculate()

        nc = netCDF(self.get_nc(), 'r')

        # read reciprocal bloch function
        npw = nc.variables['NumberPlaneWavesKpoint'][:]
        bf = nc.variables['WaveFunction'][kpt, spin, band]
        wflist = np.zeros(npw[kpt], np.complex)
        wflist.real = bf[0:npw[kpt], 1]
        wflist.imag = bf[0:npw[kpt], 0]

        nc.close()

        return wflist

    def get_reciprocal_fft_index(self, kpt=0):
        '''return the Wave Function FFT Index'''
        
        nc = netCDF(self.get_nc(), 'r')
        recind = nc.variables['WaveFunctionFFTindex'][kpt, :, :]
        nc.close()
        return recind

    def get_ensemble_coefficients(self):
        'returns exchange correlation ensemble coefficients'
        
        # adapted from ASE/dacapo.py
        # def GetEnsembleCoefficients(self):
        #     self.Calculate()
        #     E = self.GetPotentialEnergy()
        #     xc = self.GetNetCDFEntry('EnsembleXCEnergies')
        #     Exc = xc[0]
        #     exc_c = self.GetNetCDFEntry('EvaluateCorrelationEnergy')
        #     exc_e =  self.GetNetCDFEntry('EvaluateExchangeEnergy')
        #     exc = exc_c + exc_e
        #     if self.GetXCFunctional() == 'RPBE':
        #             Exc = exc[-1][-1]
        # 
        #     E0 = xc[1]       # Fx = 0
        # 
        #     diff0 = xc[2] # - Exc
        #     diff1 = xc[3] # - Exc
        #     diff2 = xc[4] # - Exc
        #     coefs = (E + E0 - Exc,diff0-E0 ,diff1-E0,diff2-E0)
        #     print 'ensemble: (%.9f, %.9f, %.9f, %.9f)'% coefs
        #     return num.array(coefs)
        if self.calculation_required():
            self.calculate()

        E = self.get_potential_energy()
        nc = netCDF(self.get_nc(), 'r')
        if 'EnsembleXCEnergies' in nc.variables:
            v = nc.variables['EnsembleXCEnergies']
            xc = v[:]

        EXC = xc[0]

        if 'EvaluateCorrelationEnergy' in nc.variables:
            v = nc.variables['EvaluateCorrelationEnergy']
            exc_c = v[:]
            
        if 'EvaluateExchangeEnergy' in nc.variables:
            v = nc.variables['EvaluateExchangeEnergy']
            exc_e = v[:]

        exc = exc_c + exc_e

        if self.get_xc == 'RPBE':
            EXC = exc[-1][-1]
            
        E0 = xc[1]    # Fx = 0

        diff0 = xc[2] # - Exc
        diff1 = xc[3] # - Exc
        diff2 = xc[4] # - Exc
        coefs = (E + E0 - EXC, diff0-E0, diff1-E0, diff2-E0)
        log.info('ensemble: (%.9f, %.9f, %.9f, %.9f)'% coefs)
        return np.array(coefs)

    def get_pseudo_wave_function(self, band=0, kpt=0, spin=0, pad=True):

        '''return the pseudo wavefunction'''
        
        # pad=True does nothing here.
        if self.calculation_required():
            self.calculate()

        ibz = self.get_ibz_kpoints()

        #get the reciprocal bloch function
        wflist = self.get_reciprocal_bloch_function(band=band,
                                                    kpt=kpt,
                                                    spin=spin)
        # wflist == Reciprocal Bloch Function
 
        recind = self. get_reciprocal_fft_index(kpt)
        grids = self.get_fftgrid()
        softgrid = grids['soft']
        
        # GetReciprocalBlochFunctionGrid
        wfrec = np.zeros((softgrid), np.complex) 

        for i in xrange(len(wflist)):
            wfrec[recind[0, i]-1,
                  recind[1, i]-1,
                  recind[2, i]-1] = wflist[i]

        # calculate Bloch Function
        wf = wfrec.copy() 
        dim = wf.shape
        for i in range(len(dim)):
            wf = np.fft.fft(wf, dim[i], axis=i)

        #now the phase function to get the bloch phase
        basis = self.get_atoms().get_cell()
        kpoint = np.dot(ibz[kpt], basis) #coordinates of relevant
                                         #kpoint in cartesian
                                         #coordinates
        def phasefunction(coor):
            'return phasefunction'
            pf = np.exp(1.0j*np.dot(kpoint, coor))
            return pf

        # Calculating the Bloch phase at the origin (0,0,0) of the grid
        origin = np.array([0., 0., 0.])
        blochphase = phasefunction(origin)
        spatialshape = wf.shape[-len(basis):]
        gridunitvectors = np.array(map(lambda unitvector,
                                       shape:unitvector/shape,
                                       basis,
                                       spatialshape))

        for dim in range(len(spatialshape)):
            # Multiplying with the phase at the origin
            deltaphase = phasefunction(gridunitvectors[dim])
            # and calculating phase difference between each point
            newphase = np.fromfunction(lambda i, phase=deltaphase:phase**i,
                                     (spatialshape[dim],))
            blochphase = np.multiply.outer(blochphase, newphase)

        return blochphase*wf

    def get_wave_function(self, band=0, kpt=0, spin=0):
        '''return the wave function. This is the pseudo wave function
        divided by volume.'''
        
        pwf = self.get_pseudo_wave_function(band=band,
                                            kpt=kpt,
                                            spin=spin,
                                            pad=True)
        vol = self.get_atoms().get_volume()
        fftgrids = self.get_fftgrid()
        softgrid = fftgrids['soft']
    
        x, y, z = self.get_ucgrid((softgrid))

        return x, y, z, pwf/np.sqrt(vol)

    def strip(self):
        '''remove all large memory nc variables not needed for
        anything I use very often. 
        '''
        self.delete_ncattdimvar(self.nc,
                                ncdims=['max_projectors_per_atom'],
                                ncvars=['WaveFunction',
                                        'WaveFunctionFFTindex',
                                        'NumberOfNLProjectors',
                                        'NLProjectorPsi',
                                        'TypeNLProjector1',
                                        'NumberofNLProjectors',
                                        'PartialCoreDensity',
                                        'ChargeDensity',
                                        'ElectrostaticPotential',
                                        'StructureFactor'])

# shortcut function names
Jacapo.get_cd = Jacapo.get_charge_density
Jacapo.get_wf = Jacapo.get_wave_function
Jacapo.get_esp = Jacapo.get_electrostatic_potential
Jacapo.get_occ = Jacapo.get_occupation_numbers
Jacapo.get_ef = Jacapo.get_fermi_level
Jacapo.get_number_of_bands = Jacapo.get_nbands
Jacapo.get_electronic_temperature = Jacapo.get_ft
Jacapo.get_number_of_electrons = Jacapo.get_valence

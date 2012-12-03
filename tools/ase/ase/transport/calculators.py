import numpy as np

from numpy import linalg
from ase.transport.selfenergy import LeadSelfEnergy, BoxProbe
from ase.transport.greenfunction import GreenFunction
from ase.transport.tools import subdiagonalize, cutcoupling, tri2full, dagger,\
    rotate_matrix


class TransportCalculator:
    """Determine transport properties of a device sandwiched between
    two semi-infinite leads using a Green function method.
    """

    def __init__(self, **kwargs):
        """Create the transport calculator.

        Parameters
        ==========
        h : (N, N) ndarray
            Hamiltonian matrix for the central region. 
        s : {None, (N, N) ndarray}, optional
            Overlap matrix for the central region. 
            Use None for an orthonormal basis.
        h1 : (N1, N1) ndarray
            Hamiltonian matrix for lead1.
        h2 : {None, (N2, N2) ndarray}, optional
            Hamiltonian matrix for lead2. You may use None if lead1 and lead2 
            are identical.
        s1 : {None, (N1, N1) ndarray}, optional
            Overlap matrix for lead1. Use None for an orthonomormal basis.
        hc1 : {None, (N1, N) ndarray}, optional
            Hamiltonian coupling matrix between the first principal
            layer in lead1 and the central region.
        hc2 : {None, (N2, N} ndarray), optional
            Hamiltonian coupling matrix between the first principal
            layer in lead2 and the central region.
        sc1 : {None, (N1, N) ndarray}, optional  
            Overlap coupling matrix between the first principal
            layer in lead1 and the central region.
        sc2 : {None, (N2, N) ndarray}, optional  
            Overlap coupling matrix between the first principal
            layer in lead2 and the central region.
        energies : {None, array_like}, optional
            Energy points for which calculated transport properties are
            evaluated.
        eta : {1.0e-5, float}, optional
            Infinitesimal for the central region Green function. 
        eta1/eta2 : {1.0e-5, float}, optional
            Infinitesimal for lead1/lead2 Green function.
        align_bf : {None, int}, optional
            Use align_bf=m to shift the central region 
            by a constant potential such that the m'th onsite element
            in the central region is aligned to the m'th onsite element
            in lead1 principal layer.
        logfile : {None, str}, optional 
            Write a logfile to file with name `logfile`.
            Use '-' to write to std out.
        eigenchannels: {0, int}, optional
            Number of eigenchannel transmission coefficients to 
            calculate. 
        pdos : {None, (N,) array_like}, optional
            Specify which basis functions to calculate the
            projected density of states for.
        dos : {False, bool}, optional
            The total density of states of the central region.
        box: XXX
            YYY
            
        If hc1/hc2 are None, they are assumed to be identical to
        the coupling matrix elements between neareste neighbor 
        principal layers in lead1/lead2.

        Examples
        ========
        >>> import numpy as np
        >>> h = np.array((0,)).reshape((1,1))
        >>> h1 = np.array((0, -1, -1, 0)).reshape(2,2)
        >>> energies = np.arange(-3, 3, 0.1)
        >>> calc = TransportCalculator(h=h, h1=h1, energies=energies)
        >>> T = calc.get_transmission()

        """
        
        # The default values for all extra keywords
        self.input_parameters = {'energies': None,
                                 'h': None,
                                 'h1': None,
                                 'h2': None,
                                 's': None,
                                 's1': None,
                                 's2': None,
                                 'hc1': None,
                                 'hc2': None,
                                 'sc1': None,
                                 'sc2': None,
                                 'box': None,
                                 'align_bf': None,
                                 'eta1': 1e-5,
                                 'eta2': 1e-5,
                                 'eta': 1e-5,
                                 'logfile': None, # '-',
                                 'eigenchannels': 0,
                                 'dos': False,
                                 'pdos': [],
                                 }
        self.initialized = False # Changed Hamiltonians?
        self.uptodate = False # Changed energy grid?
        self.set(**kwargs)

    def set(self, **kwargs):
        for key in kwargs:
            if key in ['h', 'h1', 'h2', 'hc1', 'hc2',
                       's', 's1', 's2', 'sc1', 'sc2',
                       'eta', 'eta1', 'eta2', 'align_bf', 'box']:
                self.initialized = False
                self.uptodate = False
                break
            elif key in ['energies', 'eigenchannels', 'dos', 'pdos']:
                self.uptodate = False
            elif key not in self.input_parameters:
                raise KeyError, '\'%s\' not a vaild keyword' % key

        self.input_parameters.update(kwargs)
        log = self.input_parameters['logfile']
        if log is None:
            class Trash:
                def write(self, s):
                    pass
                def flush(self):
                    pass
            self.log = Trash()
        elif log == '-':
            from sys import stdout
            self.log = stdout
        elif 'logfile' in kwargs:
            self.log = open(log, 'w')

    def initialize(self):
        if self.initialized:
            return

        print >> self.log, '# Initializing calculator...'

        p = self.input_parameters
        if p['s'] == None:
            p['s'] = np.identity(len(p['h']))
        
        identical_leads = False
        if p['h2'] == None:   
            p['h2'] = p['h1'] # Lead2 is idendical to lead1
            identical_leads = True
 
        if p['s1'] == None: 
            p['s1'] = np.identity(len(p['h1']))
       
        if p['s2'] == None and not identical_leads:
            p['s2'] = np.identity(len(p['h2'])) # Orthonormal basis for lead 2
        else: # Lead2 is idendical to lead1
            p['s2'] = p['s1']

           
        h_mm = p['h']
        s_mm = p['s']
        pl1 = len(p['h1']) / 2
        pl2 = len(p['h2']) / 2
        h1_ii = p['h1'][:pl1, :pl1]
        h1_ij = p['h1'][:pl1, pl1:2 * pl1]
        s1_ii = p['s1'][:pl1, :pl1]
        s1_ij = p['s1'][:pl1, pl1:2 * pl1]
        h2_ii = p['h2'][:pl2, :pl2]
        h2_ij = p['h2'][pl2: 2 * pl2, :pl2]
        s2_ii = p['s2'][:pl2, :pl2]
        s2_ij = p['s2'][pl2: 2 * pl2, :pl2]
        
        if p['hc1'] is None:
            nbf = len(h_mm)
            h1_im = np.zeros((pl1, nbf), complex)
            s1_im = np.zeros((pl1, nbf), complex)
            h1_im[:pl1, :pl1] = h1_ij
            s1_im[:pl1, :pl1] = s1_ij
            p['hc1'] = h1_im
            p['sc1'] = s1_im
        else:
            h1_im = p['hc1']
            if p['sc1'] is not None:
                s1_im = p['sc1']
            else:
                s1_im = np.zeros(h1_im.shape, complex)
                p['sc1'] = s1_im

        if p['hc2'] is None:
            h2_im = np.zeros((pl2, nbf), complex)
            s2_im = np.zeros((pl2, nbf), complex)
            h2_im[-pl2:, -pl2:] = h2_ij
            s2_im[-pl2:, -pl2:] = s2_ij
            p['hc2'] = h2_im
            p['sc2'] = s2_im
        else:
            h2_im = p['hc2']
            if p['sc2'] is not None:
                s2_im = p['sc2']
            else:
                s2_im = np.zeros(h2_im.shape, complex)
                p['sc2'] = s2_im

        align_bf = p['align_bf']
        if align_bf != None:
            diff = (h_mm[align_bf, align_bf] - h1_ii[align_bf, align_bf]) \
                   / s_mm[align_bf, align_bf]
            print >> self.log, '# Aligning scat. H to left lead H. diff=', diff
            h_mm -= diff * s_mm

        # Setup lead self-energies
        # All infinitesimals must be > 0 
        assert np.all(np.array((p['eta'], p['eta1'], p['eta2'])) > 0.0)
        self.selfenergies = [LeadSelfEnergy((h1_ii, s1_ii), 
                                            (h1_ij, s1_ij),
                                            (h1_im, s1_im),
                                            p['eta1']),
                             LeadSelfEnergy((h2_ii, s2_ii), 
                                            (h2_ij, s2_ij),
                                            (h2_im, s2_im),
                                            p['eta2'])]
        box = p['box']
        if box is not None:
            print 'Using box probe!'
            self.selfenergies.append(
                BoxProbe(eta=box[0], a=box[1], b=box[2], energies=box[3],
                         S=s_mm, T=0.3))
        
        #setup scattering green function
        self.greenfunction = GreenFunction(selfenergies=self.selfenergies,
                                           H=h_mm,
                                           S=s_mm,
                                           eta=p['eta'])

        self.initialized = True
    
    def update(self):
        if self.uptodate:
            return
        
        p = self.input_parameters
        self.energies = p['energies']
        nepts = len(self.energies)
        nchan = p['eigenchannels']
        pdos = p['pdos']
        self.T_e = np.empty(nepts)
        if p['dos']:
            self.dos_e = np.empty(nepts)
        if pdos != []:
            self.pdos_ne = np.empty((len(pdos), nepts))
        if nchan > 0:
            self.eigenchannels_ne = np.empty((nchan, nepts))

        for e, energy in enumerate(self.energies):
            Ginv_mm = self.greenfunction.retarded(energy, inverse=True)
            lambda1_mm = self.selfenergies[0].get_lambda(energy)
            lambda2_mm = self.selfenergies[1].get_lambda(energy)
            a_mm = linalg.solve(Ginv_mm, lambda1_mm)
            b_mm = linalg.solve(dagger(Ginv_mm), lambda2_mm)
            T_mm = np.dot(a_mm, b_mm)
            if nchan > 0:
                t_n = linalg.eigvals(T_mm).real
                self.eigenchannels_ne[:, e] = np.sort(t_n)[-nchan:]
                self.T_e[e] = np.sum(t_n)
            else:
                self.T_e[e] = np.trace(T_mm).real

            print >> self.log, energy, self.T_e[e]
            self.log.flush()

            if p['dos']:
                self.dos_e[e] = self.greenfunction.dos(energy)

            if pdos != []:
                self.pdos_ne[:, e] = np.take(self.greenfunction.pdos(energy),
                                             pdos)
        
        self.uptodate = True

    def print_pl_convergence(self):
        self.initialize()
        pl1 = len(self.input_parameters['h1']) / 2
        
        h_ii = self.selfenergies[0].h_ii
        s_ii = self.selfenergies[0].s_ii
        ha_ii = self.greenfunction.H[:pl1, :pl1]
        sa_ii = self.greenfunction.S[:pl1, :pl1]
        c1 = np.abs(h_ii - ha_ii).max()
        c2 = np.abs(s_ii - sa_ii).max()
        print 'Conv (h,s)=%.2e, %2.e' % (c1, c2)

    def plot_pl_convergence(self):
        self.initialize()
        pl1 = len(self.input_parameters['h1']) / 2       
        hlead = self.selfenergies[0].h_ii.real.diagonal()
        hprincipal = self.greenfunction.H.real.diagonal[:pl1]

        import pylab as pl
        pl.plot(hlead, label='lead')
        pl.plot(hprincipal, label='principal layer')
        pl.axis('tight')
        pl.show()

    def get_transmission(self):
        self.initialize()
        self.update()
        return self.T_e

    def get_dos(self):
        self.initialize()
        self.update()
        return self.dos_e

    def get_eigenchannels(self, n=None):
        """Get ``n`` first eigenchannels."""
        self.initialize()
        self.update()
        if n is None:
            n = self.input_parameters['eigenchannels']
        return self.eigenchannels_ne[:n]

    def get_pdos(self):
        self.initialize()
        self.update()
        return self.pdos_ne

    def subdiagonalize_bfs(self, bfs, apply=False):
        self.initialize()
        bfs = np.array(bfs)
        p = self.input_parameters
        h_mm = p['h']
        s_mm = p['s']
        ht_mm, st_mm, c_mm, e_m = subdiagonalize(h_mm, s_mm, bfs)
        if apply:
            self.uptodate = False
            h_mm[:] = ht_mm 
            s_mm[:] = st_mm 
            # Rotate coupling between lead and central region
            for alpha, sigma in enumerate(self.selfenergies):
                sigma.h_im[:] = np.dot(sigma.h_im, c_mm)
                sigma.s_im[:] = np.dot(sigma.s_im, c_mm)
        
        c_mm = np.take(c_mm, bfs, axis=0)
        c_mm = np.take(c_mm, bfs, axis=1)
        return ht_mm, st_mm, e_m, c_mm

    def cutcoupling_bfs(self, bfs, apply=False):
        self.initialize()
        bfs = np.array(bfs)
        p = self.input_parameters
        h_pp = p['h'].copy()
        s_pp = p['s'].copy()
        cutcoupling(h_pp, s_pp, bfs)
        if apply:
            self.uptodate = False
            p['h'][:] = h_pp
            p['s'][:] = s_pp
            for alpha, sigma in enumerate(self.selfenergies):
                for m in bfs:
                    sigma.h_im[:, m] = 0.0
                    sigma.s_im[:, m] = 0.0
        return h_pp, s_pp

    def lowdin_rotation(self, apply=False):
        p = self.input_parameters
        h_mm = p['h']
        s_mm = p['s']
        eig, rot_mm = linalg.eigh(s_mm)
        eig = np.abs(eig)
        rot_mm = np.dot(rot_mm / np.sqrt(eig), dagger(rot_mm))
        if apply:
            self.uptodate = False
            h_mm[:] = rotate_matrix(h_mm, rot_mm) # rotate C region
            s_mm[:] = rotate_matrix(s_mm, rot_mm)
            for alpha, sigma in enumerate(self.selfenergies):
                sigma.h_im[:] = np.dot(sigma.h_im, rot_mm) # rotate L-C coupl.
                sigma.s_im[:] = np.dot(sigma.s_im, rot_mm)

        return rot_mm

    def get_left_channels(self, energy, nchan=1):
        self.initialize()
        g_s_ii = self.greenfunction.retarded(energy)
        lambda_l_ii = self.selfenergies[0].get_lambda(energy)
        lambda_r_ii = self.selfenergies[1].get_lambda(energy)

        if self.greenfunction.S is None:
            s_s_qsrt_ii = s_s_isqrt = np.identity(len(g_s_ii))
        else:
            s_mm = self.greenfunction.S
            s_s_i, s_s_ii = linalg.eig(s_mm)
            s_s_i = np.abs(s_s_i)
            s_s_sqrt_i = np.sqrt(s_s_i) # sqrt of eigenvalues  
            s_s_sqrt_ii = np.dot(s_s_ii * s_s_sqrt_i, dagger(s_s_ii))
            s_s_isqrt_ii = np.dot(s_s_ii / s_s_sqrt_i, dagger(s_s_ii))

        lambdab_r_ii = np.dot(np.dot(s_s_isqrt_ii, lambda_r_ii),s_s_isqrt_ii)
        a_l_ii = np.dot(np.dot(g_s_ii, lambda_l_ii), dagger(g_s_ii))
        ab_l_ii = np.dot(np.dot(s_s_sqrt_ii, a_l_ii), s_s_sqrt_ii)
        lambda_i, u_ii = linalg.eig(ab_l_ii)
        ut_ii = np.sqrt(lambda_i / (2.0 * np.pi)) * u_ii
        m_ii = 2 * np.pi * np.dot(np.dot(dagger(ut_ii), lambdab_r_ii),ut_ii)
        T_i,c_in = linalg.eig(m_ii)
        T_i = np.abs(T_i)
        
        channels = np.argsort(-T_i)[:nchan]
        c_in = np.take(c_in, channels, axis=1)
        T_n = np.take(T_i, channels)
        v_in = np.dot(np.dot(s_s_isqrt_ii, ut_ii), c_in)

        return T_n, v_in

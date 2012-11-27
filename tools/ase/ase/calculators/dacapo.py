import numpy as np

from ase.old import OldASEListOfAtomsWrapper

try:
    import Numeric as num
except ImportError:
    pass

def np2num(a, typecode=None):
    if num.__version__ > '23.8':
        return num.array(a, typecode)
    if typecode is None:
        typecode = num.Float
    b = num.fromstring(a.tostring(), typecode)
    b.shape = a.shape
    return b

def restart(filename, **kwargs):
    calc = Dacapo(filename, **kwargs)
    atoms = calc.get_atoms()
    return atoms, calc

class Dacapo:
    def __init__(self, filename=None, stay_alive=False, stress=False,
                 **kwargs):

        self.kwargs = kwargs
        self.stay_alive = stay_alive
        self.stress = stress
        
        if filename is not None:
            from Dacapo import Dacapo
            self.loa = Dacapo.ReadAtoms(filename, **kwargs)
            self.calc = self.loa.GetCalculator()
        else:
            self.loa = None
            self.calc = None

        self.pps = []
        
    def set_pp(self, Z, path):
        self.pps.append((Z, path))

    def set_txt(self, txt):
        if self.calc is None:
            self.kwargs['txtout'] = txt
        else:
            self.calc.SetTxtFile(txt)

    def set_nc(self, nc):
        if self.calc is None:
            self.kwargs['out'] = nc
        else:
            self.calc.SetNetCDFFile(nc)

    def update(self, atoms):
        from Dacapo import Dacapo
        if self.calc is None:
            if 'nbands' not in self.kwargs:
                n = sum([valence[atom.symbol] for atom in atoms])
                self.kwargs['nbands'] = int(n * 0.65) + 4

            magmoms = atoms.get_initial_magnetic_moments()
            if magmoms.any():
                self.kwargs['spinpol'] = True

            self.calc = Dacapo(**self.kwargs)

            if self.stay_alive:
                self.calc.StayAliveOn()
            else:
                self.calc.StayAliveOff()

            if self.stress:
                self.calc.CalculateStress()

            for Z, path in self.pps:
                self.calc.SetPseudoPotential(Z, path)

        if self.loa is None:
            from ASE import Atom, ListOfAtoms
            numbers = atoms.get_atomic_numbers()
            positions = atoms.get_positions()
            magmoms = atoms.get_initial_magnetic_moments()
            self.loa = ListOfAtoms([Atom(Z=numbers[a],
                                         position=positions[a],
                                         magmom=magmoms[a])
                                    for a in range(len(atoms))],
                                   cell=np2num(atoms.get_cell()),
                                   periodic=tuple(atoms.get_pbc()))
            self.loa.SetCalculator(self.calc)
        else:
            self.loa.SetCartesianPositions(np2num(atoms.get_positions()))
            self.loa.SetUnitCell(np2num(atoms.get_cell()), fix=True)
            
    def get_atoms(self):
        atoms = OldASEListOfAtomsWrapper(self.loa).copy()
        atoms.set_calculator(self)
        return atoms
    
    def get_potential_energy(self, atoms):
        self.update(atoms)
        return self.calc.GetPotentialEnergy()

    def get_forces(self, atoms):
        self.update(atoms)
        return np.array(self.calc.GetCartesianForces())

    def get_stress(self, atoms):
        self.update(atoms)
        stress = np.array(self.calc.GetStress())
        if stress.ndim == 2:
            return stress.ravel()[[0, 4, 8, 5, 2, 1]]
        else:
            return stress

    def calculation_required(self, atoms, quantities):
        if self.calc is None:
            return True

        if atoms != self.get_atoms():
            return True

        return False
        
    def get_number_of_bands(self):
        return self.calc.GetNumberOfBands()

    def get_k_point_weights(self):
        return np.array(self.calc.GetIBZKPointWeights())

    def get_number_of_spins(self):
        return 1 + int(self.calc.GetSpinPolarized())

    def get_eigenvalues(self, kpt=0, spin=0):
        return np.array(self.calc.GetEigenvalues(kpt, spin))

    def get_fermi_level(self):
        return self.calc.GetFermiLevel()

    def get_magnetic_moment(self):
        return self.calc.GetMagneticMoment()

    def get_number_of_electrons(self):
        return self.calc.GetValenceElectrons()

    def get_number_of_grid_points(self):
        return np.array(self.get_pseudo_wave_function(0, 0, 0).shape)

    def get_pseudo_density(self, spin=0):
        return np.array(self.calc.GetDensityArray(s))
    
    def get_pseudo_wave_function(self, band=0, kpt=0, spin=0, pad=True):
        kpt_c = self.get_bz_k_points()[kpt]
        state = self.calc.GetElectronicStates().GetState(band=band, spin=spin,
                                                         kptindex=kpt)

        # Get wf, without bloch phase (Phase = True doesn't do anything!)
        wave = state.GetWavefunctionOnGrid(phase=False)

        # Add bloch phase if this is not the Gamma point
        if np.all(kpt_c == 0):
            return wave
        coord = state.GetCoordinates()
        phase = coord[0] * kpt_c[0] + coord[1] * kpt_c[1] + coord[2] * kpt_c[2]
        return np.array(wave) * np.exp(-2.j * np.pi * phase) # sign! XXX

        #return np.array(self.calc.GetWaveFunctionArray(n, k, s)) # No phase!

    def get_bz_k_points(self):
        return np.array(self.calc.GetBZKPoints())

    def get_ibz_k_points(self):
        return np.array(self.calc.GetIBZKPoints())

    def get_wannier_localization_matrix(self, nbands, dirG, kpoint,
                                        nextkpoint, G_I, spin):
        return np.array(self.calc.GetWannierLocalizationMatrix(
            G_I=G_I.tolist(), nbands=nbands, dirG=dirG.tolist(),
            kpoint=kpoint, nextkpoint=nextkpoint, spin=spin))
    
    def initial_wannier(self, initialwannier, kpointgrid, fixedstates,
                        edf, spin):
        # Use initial guess to determine U and C
        init = self.calc.InitialWannier(initialwannier, self.atoms,
                                        np2num(kpointgrid, num.Int))

        states = self.calc.GetElectronicStates()
        waves = [[state.GetWaveFunction()
                  for state in states.GetStatesKPoint(k, spin)]
                 for k in self.calc.GetIBZKPoints()] 

        init.SetupMMatrix(waves, self.calc.GetBZKPoints())
        c, U = init.GetListOfCoefficientsAndRotationMatrices(
            (self.calc.GetNumberOfBands(), fixedstates, edf))
        U = np.array(U)
        for k in range(len(c)):
            c[k] = np.array(c[k])
        return c, U

valence = {
'H':   1,
'B':   3,
'C':   4,
'N':   5,
'O':   6,
'Li':  1,
'Na':  1,
'K':   9,
'Mg':  8,
'Ca': 10,
'Sr': 10,
'Al':  3,
'Ga': 13,
'Sc': 11,
'Ti': 12,
'V':  13,
'Cr': 14,
'Mn':  7,
'Fe':  8,
'Co':  9,
'Ni': 10,
'Cu': 11,
'Zn': 12,
'Y':  11,
'Zr': 12,
'Nb': 13,
'Mo':  6,
'Ru':  8,
'Rh':  9,
'Pd': 10,
'Ag': 11,
'Cd': 12,
'Ir': 9,
'Pt': 10,
'Au': 11,
}

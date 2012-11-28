import commands, os, string, tempfile, shutil
from ase import write
from ase.units import Bohr

class Bader:
    '''class for running bader analysis and extracting data from it.

    The class runs bader, extracts the charge density and outputs it
    to a cube file. Then you call different functions of the class to
    extract the charges, volumes, etc...

    ACF.dat contains the coordinates of each atom, the charge
    associated with it according to Bader partitioning, percentage of
    the whole according to Bader partitioning and the minimum distance
    to the surface. This distance should be compared to maximum
    cut-off radius for the core region if pseudo potentials have been
    used.

    BCF.dat contains the coordinates of each Bader maxima, the charge
    within that volume, the nearest atom and the distance to that
    atom.

    AtomVolumes.dat contains the number of each volume that has been
    assigned to each atom. These numbers correspond to the number of
    the BvAtxxxx.dat files.

    The options for the executable are::

        bader [ -c bader | voronoi ]
              [ -n bader | voronoi ]
              [ -b neargrid | ongrid ]
              [ -r refine_edge_iterations ]
              [ -ref reference_charge ]
              [ -p all_atom | all_bader ]
              [ -p sel_atom | sel_bader ] [volume list]
              [ -p atom_index | bader_index ]
              [ -i cube | chgcar ]
              [ -h ] [ -v ]
              chargefile

    References:
    
    G. Henkelman, A. Arnaldsson, and H. Jonsson, A fast and robust
    algorithm for Bader decomposition of charge density,
    Comput. Mater. Sci. 36 254-360 (2006).

    E. Sanville, S. D. Kenny, R. Smith, and G. Henkelman An improved
    grid-based algorithm for Bader charge allocation,
    J. Comp. Chem. 28 899-908 (2007).

    W. Tang, E. Sanville, and G. Henkelman A grid-based Bader analysis
    algorithm without lattice bias, J. Phys.: Condens. Matter 21
    084204 (2009).
    '''
    def __init__(self,atoms):
        '''
 
        '''
        self.atoms = atoms

        #get density and write cube file
        calc = atoms.get_calculator()
        ncfile = calc.get_nc()
        base,ext = os.path.splitext(ncfile)

        x,y,z,density = calc.get_charge_density()
        cubefile = base + '_charge_density.cube'
        self.densityfile = cubefile

        if not os.path.exists(cubefile):
            write(cubefile, atoms,data=density*Bohr**3)
        
        #cmd to run for bader analysis. check if output exists so we
        #don't run this too often.
        acf_file = base + '_ACF.dat'
        if not os.path.exists(acf_file):
            #mk tempdir
            tempdir = tempfile.mkdtemp()
            
            cwd = os.getcwd()
            os.chdir(tempdir)
            
            cmd = 'bader %s' % abscubefile
            status,output = commands.getstatusoutput(cmd)
            
            if status != 0:
                print output

            shutil.copy2('ACF.dat',os.path.join(cwd,acf_file))
            
            os.chdir(cwd)
            shutil.rmtree(tempdir)

        self.charges = []
        self.volumes = []

        #now parse the output
        f = open(acf_file,'r')
        #skip 2 lines
        f.readline()
        f.readline()

        for i,atom in enumerate(self.atoms):
            line = f.readline()
            fields = line.split()
            n = int(fields[0])
            x = float(fields[1])
            y = float(fields[2])
            z = float(fields[3])
            chg = float(fields[4])
            mindist = float(fields[5])
            vol = float(fields[6])

            self.charges.append(chg)
            self.volumes.append(vol)

        f.close()

    def get_bader_charges(self):
        return self.charges

    def get_bader_volumes(self):
        'return volumes in Ang**3'
        return [x*Bohr**3 for x in self.volumes]

    def write_atom_volume(self,atomlist):
        '''write bader atom volumes to cube files.
        atomlist = [0,2] #for example

        -p sel_atom Write the selected atomic volumes, read from the
        subsequent list of volumes.
        '''
        alist = string.join([str(x) for x in atomlist],' ')
        cmd = 'bader -p sel_atom %s %s' % (alist,self.densityfile)
        print cmd
        os.system(cmd)
        
    def write_bader_volume(self,atomlist):
        """write bader atom volumes to cube files.

        ::
        
          atomlist = [0,2] #  for example
          
        -p sel_bader Write the selected Bader volumes, read from the
        subsequent list of volumes.
        """
        alist = string.join([str(x) for x in atomlist],' ')
        cmd = 'bader -p sel_bader %s %s' % (alist,self.densityfile)
        print cmd
        os.system(cmd)

    def write_atom_index(self):
        ''' -p atom_index Write the atomic volume index to a charge
        density file.
        '''
        cmd = 'bader -p atom_index %s' % (self.densityfile)
        print cmd
        os.system(cmd)

    def write_bader_index(self):
        '''
        -p bader_index Write the Bader volume index to a charge
        density file.
        '''
        cmd = 'bader -p bader_index %s' % (self.densityfile)
        print cmd
        os.system(cmd)

    def write_all_atom(self):
        '''
        -p all_atom Combine all volumes associated with an atom and
        write to file. This is done for all atoms and written to files
        named BvAtxxxx.dat. The volumes associated with atoms are
        those for which the maximum in charge density within the
        volume is closest to the atom.
        '''
        cmd = 'bader -p all_atom %s' % (self.densityfile)
        print cmd
        os.system(cmd)

    def write_all_bader(self):
        '''
        -p all_bader Write all Bader volumes (containing charge above
        threshold of 0.0001) to a file. The charge distribution in
        each volume is written to a separate file, named
        Bvolxxxx.dat. It will either be of a CHGCAR format or a CUBE
        file format, depending on the format of the initial charge
        density file. These files can be quite large, so this option
        should be used with caution.
        '''
        cmd = 'bader -p all_bader %s' % (self.densityfile)
        print cmd
        os.system(cmd)
        
if __name__ == '__main__':

    from Jacapo import *

    atoms = Jacapo.read_atoms('ethylene.nc')

    b = Bader(atoms)

    print b.get_bader_charges()
    print b.get_bader_volumes()
    b.write_atom_volume([3,4])

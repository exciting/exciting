from ase.atoms import Atoms
from ase.quaternions import Quaternions
from ase.parallel import paropen

def read_lammps_dump(fileobj, index=-1):
    """Method which reads a LAMMPS dump file."""
    if isinstance(fileobj, str):
        f = paropen(fileobj)
    else:
        f = fileobj

    # load everything into memory
    lines = f.readlines()

    natoms = 0
    images = []

    while len(lines) > natoms:
        line = lines.pop(0)

        if 'ITEM: TIMESTEP' in line:
            n_atoms = 0
            lo = [] ; hi = [] ; tilt = []
            id = [] ; type = []
            positions = []
            velocities = [] 
            forces = []
            quaternions = []

        if 'ITEM: NUMBER OF ATOMS' in line:
            line = lines.pop(0)
            natoms = int(line.split()[0])
            
        if 'ITEM: BOX BOUNDS' in line:
            # save labels behind "ITEM: BOX BOUNDS" in triclinic case (>=lammps-7Jul09)
            tilt_items = line.split()[3:]
            for i in range(3):
                line = lines.pop(0)
                fields = line.split()
                lo.append(float(fields[0]))
                hi.append(float(fields[1]))
                if (len(fields) >= 3):
                    tilt.append(float(fields[2]))

            # determine cell tilt (triclinic case!)
            if (len(tilt) >= 3):
                # for >=lammps-7Jul09 use labels behind "ITEM: BOX BOUNDS" to assign tilt (vector) elements ...
                if (len(tilt_items) >= 3):
                    xy = tilt[tilt_items.index('xy')]
                    xz = tilt[tilt_items.index('xz')]
                    yz = tilt[tilt_items.index('yz')]
                # ... otherwise assume default order in 3rd column (if the latter was present)
                else:
                    xy = tilt[0]
                    xz = tilt[1]
                    yz = tilt[2]
            else:
                xy = xz = yz = 0
            xhilo = (hi[0] - lo[0]) - xy - xz
            yhilo = (hi[1] - lo[1]) - yz
            zhilo = (hi[2] - lo[2])

            cell = [[xhilo,0,0],[xy,yhilo,0],[xz,yz,zhilo]]

        def add_quantity(fields, var, labels):
            for label in labels:
                if label not in atom_attributes:
                    return
            var.append([float(fields[atom_attributes[label]])
                        for label in labels])
                
        if 'ITEM: ATOMS' in line:
            # (reliably) identify values by labels behind "ITEM: ATOMS" - requires >=lammps-7Jul09
            # create corresponding index dictionary before iterating over atoms to (hopefully) speed up lookups...
            atom_attributes = {}
            for (i, x) in enumerate(line.split()[2:]):
                atom_attributes[x] = i
            for n in range(natoms):
                line = lines.pop(0)
                fields = line.split()
                id.append( int(fields[atom_attributes['id']]) )
                type.append( int(fields[atom_attributes['type']]) )
                add_quantity(fields, positions, ['x', 'y', 'z'])
                add_quantity(fields, velocities, ['vx', 'vy', 'vz'])
                add_quantity(fields, forces, ['fx', 'fy', 'fz'])
                add_quantity(fields, quaternions, ['c_q[1]', 'c_q[2]',
                                                   'c_q[3]', 'c_q[4]'])

            if len(quaternions):
                images.append(Quaternions(symbols=type,
                                          positions=positions,
                                          cell=cell,
                                          quaternions=quaternions))
            else:
                images.append(Atoms(symbols=type,
                                    positions=positions,
                                    cell=cell))

    return images[index]


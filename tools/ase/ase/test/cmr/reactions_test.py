from ase.test import NotAvailable
import warnings
# cmr calls all available methods in ase.atoms detected by the module inspect.
# Therefore also deprecated methods are called - and we choose to silence those warnings.
warnings.filterwarnings('ignore', 'ase.atoms.*deprecated',)

try:
    import cmr
except ImportError:
    raise NotAvailable('CMR is required')

from cmr.ui import DirectoryReader

from cmr.test.examples.ase_reaction_energy import ASEReactionEnergy

# see the module for the required format of reactions definition
from ase.test.cmr.reactions import reactions
from ase.test.cmr.reactions import reference

# assure that all reactions define a reaction_id
for r in reactions:
    assert r[-1][0] == 'reaction_id'

# project id: must uniquely identify the project!
project_id = 'EMT' + ' reaction energies'

# if True, then results are uploaded to the database
database = False

# create assisting class for project with project_id,
# that allows one to convert trajectory files into
# db-files and perform analysis
re = ASEReactionEnergy(project_id, reactions, prefix='', verbose=False)

# compounds names
compounds = re.get_compounds()

# put additional fields here:
cmr_params = {'calculator': 'EMT'}
# convert all traj files in this directory to db-files
re.create_db_files(cmr_params)

# calculate the reaction energies and write the results to db-files
# named 'reaction_id.index.db'
# Each db-file defines a group (group of all compounds belonging to
# the given reaction).

# reaction energies on initial, unoptimized geometries
cmr_params = {'geometries': 'initial'}
re.make_reaction_groups(database=False, index=0, cmr_params=cmr_params)

# print
re.print_result(database=False)

# reaction energies on final, optimized geometries
cmr_params = {'geometries': 'final'}
re.make_reaction_groups(database=False, index= -1, cmr_params=cmr_params)

# print
re.print_result(database=False)


reader = DirectoryReader('.')

# retrieve all reactions (groups) with project_id and optimized geometries from the current directory
all = reader.find(name_value_list=[('db_calculator', 'group'),
                                   ('geometries', 'final')
                                  ],
                  keyword_list=[project_id, 'reaction'])

print 'reaction_id, calc, ref, calc - ref'
# compare with the reference
for r in reactions:
    reaction_id = r[-1][1]
    res = all.get('reaction_id', reaction_id)
    if res is None:
        print "Could not find reaction_id %s in reference"%str(reaction_id)
    else:
        calc = res['reaction_energy']
        ref = reference[reaction_id]
        print reaction_id, calc, ref, calc - ref
        assert abs(calc - ref) < 1e-5



# upload the created groups to the database
if database:
    re.upload_to_database()

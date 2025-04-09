# Comparative modeling by the automodel class
#
# Demonstrates how to build multi-chain models, and symmetry restraints
#
from modeller import *              # Load standard Modeller classes
from modeller.automodel import *    # Load the automodel class
import sys

log.verbose()

class MyModel(AutoModel):
    def select_atoms(self): 
        return selection(self.residue_range('38:A', '45:A'), self.residue_range('300:A', '334:A'), self.residue_range('849:B', '856:B'), self.residue_range('1111:B', '1145:B'), self.residue_range('1660:C', '1667:C'), self.residue_range('1922:C', '1956:C'), self.residue_range('2471:D', '2478:D'), self.residue_range('2733:D', '2767:D'),  self.residue_range('3282:E', '3289:E'), self.residue_range('3544:E', '3578:E'), self.residue_range('4093:F', '4100:F'), self.residue_range('4355:F', '4389:F'))

        # The same thing from chain A (required for multi-chain models):
        # return selection(self.residue_range('1:A', '2:A'))

        # Residues 4, 6, 10:
        # return selection(self.residues['4'], self.residues['6'],
        #                  self.residues['10'])

        # All residues except 1-5:
        # return selection(self) - selection(self.residue_range('1', '5'))

env = environ()
# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']
env.io.hydrogen = env.io.hetatm = env.io.water = True

# Be sure to use 'MyModel' rather than 'automodel' here!
a = MyModel(env,
            alnfile  = 'miss-residues.ali' ,     # alignment filename
            knowns   = ('6nyf2mem_del27-45','1_37helix'),    # codes of the templates
            sequence = 'SE',
            assess_methods = (assess.DOPE,
                              #soap_protein_od.Scorer(),
                              assess.GA341))              # code of the target

a.starting_model= 1                # index of the first model
a.ending_model  = 20                # index of the last model
                                   # (determines how many models to calculate)
# Very thorough VTFM optimization:
a.library_schedule = autosched.slow
# a.max_var_iterations = 300

# Thorough MD optimization:
a.md_level = refine.fast

# Repeat the whole cycle 2 times and do not stop unless obj.func. > 1E6
#a.repeat_optimization = 2
#a.max_molpdf = 1e6


a.make()                            # do comparative modeling

ok_models = [x for x in a.outputs if x['failure'] is None]

# Rank the models by DOPE score
key = 'DOPE score'
if sys.version_info[:2] == (2,3):
    # Python 2.3's sort doesn't have a 'key' argument
    ok_models.sort(lambda a,b: cmp(a[key], b[key]))
else:
    ok_models.sort(key=lambda a: a[key])

# Get top model
m = ok_models[0]
print("Top model: %s (DOPE score %.3f)" % (m['name'], m[key]))

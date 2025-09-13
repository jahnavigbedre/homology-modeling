from modeller import *
from modeller.automodel import *

# Initialize Modeller environment
env = Environ()
env.io.atom_files_directory = ['./data/templates', '.']

class MyModel(AutoModel):
    def select_atoms(self):
        # Create and return a selection of atoms
        return Selection(self)  # Return proper Selection object instead of self.atoms

a = MyModel(env,
            alnfile='data/alignment.ali',
            knowns='1crn',   # template
            sequence='Query')  # query
a.starting_model = 1
a.ending_model   = 2  # build 2 models
a.make()

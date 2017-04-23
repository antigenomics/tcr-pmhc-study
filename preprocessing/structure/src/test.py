from util import *
from os import chdir, path, makedirs
from Bio.PDB import *

chdir(path.dirname(path.realpath(__file__)))

pdb_dir = path.abspath('../tmp/pdb/')

pdb_list = PDBList()
pdb_parser = PDBParser()

# Superimposition

pdb_file_1 = pdb_list.retrieve_pdb_file("5isz", pdir=pdb_dir)
model_1 = pdb_parser.get_structure("5isz", pdb_file_1)[0]

pdb_file_2 = pdb_list.retrieve_pdb_file("4mnq", pdir=pdb_dir)
model_2 = pdb_parser.get_structure("4mnq", pdb_file_2)[0]

residues_1 = get_residues(model_1["E"], range(90, 103))

print(unpack_residues(residues_1))

residues_2 = get_residues(model_2["E"], range(89, 101))

print(unpack_residues(residues_2))

s21 = superimpose(residues_1, residues_2, "CDR3")

print(unpack_residues(s21))

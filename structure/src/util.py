import itertools
import os
import subprocess

from Bio.PDB import PDBParser
from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile


def clear_folder(d):
    for the_file in os.listdir(d):
        file_path = os.path.join(d, the_file)
        if os.path.isfile(file_path):
            os.unlink(file_path)


# (Bio)PDB methods

def get_required_chains(pdb_group):
    return set(itertools.chain(*[pdb_group[c].tolist() for c in
                                 ['chain_tcr', 'chain_antigen', 'chain_mhc_a', 'chain_mhc_b']]))


def select_model(pdb_id, pdb_file, pdb_group):
    required_chains = get_required_chains(pdb_group)

    for i, model in enumerate(PDBParser().get_structure(pdb_id, pdb_file)):
        if required_chains.issuperset(set([chain.get_id() for chain in model])):
            return i

    return -1


def fix_pdb(pdb_id, pdb_file, pdb_group):
    chains_to_retain = get_required_chains(pdb_group)
    chains_to_remove = []
    for model in PDBParser().get_structure(pdb_id, pdb_file):
        for chain in model:
            if chain.get_id() not in chains_to_retain:
                chains_to_remove.append(chain.get_id())

    fixer = PDBFixer(filename=pdb_file)

    fixer.removeChains(chainIds=chains_to_remove)

    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.removeHeterogens(True)

    PDBFile.writeFile(fixer.topology, fixer.positions, open(pdb_file, 'w'))


_trans = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
          'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
          'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
          'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}


def get_aa_code(res):
    return _trans[res.get_resname()]


def get_seq(residues):
    return "".join(get_aa_code(res) for res in residues)


def get_residues(residues, range):
    return [x for i, x in enumerate(residues) if i in range]


def calc_distance(aa1, aa2):
    return min([abs(atom1 - atom2) for atom1 in aa1 for atom2 in aa2])


def calc_distances(tcr_residues, ag_residues, tcr_range, ag_range):
    return [{'aa_tcr': get_aa_code(aa_tcr),
             'aa_ag': get_aa_code(aa_ag),
             'len_tcr': len(tcr_range),
             'len_ag': len(ag_range),
             'pos_tcr': tcr_range[i],
             'pos_ag': ag_range[j],
             'dist': calc_distance(aa_tcr, aa_ag)}
            for i, aa_tcr in enumerate(tcr_residues)
            for j, aa_ag in enumerate(ag_residues)]


# GROMACS methods

def prepare_gmx(pdb_file, gmx_prefix):
    pdb_file = os.path.abspath(pdb_file)
    gmx_prefix = os.path.abspath(gmx_prefix)

    if os.name == 'nt':  # should run via cygwin on windows
        pdb_file = pdb_file.replace('\\', '\\\\\\\\')
        gmx_prefix = gmx_prefix.replace('\\', '\\\\\\\\')

    subprocess.check_call(
            'gmx pdb2gmx -f "{}" -o "{}.gro" -p "{}.top" -n "{}" -water spce -ff oplsaa'.
                format(pdb_file, gmx_prefix, gmx_prefix, gmx_prefix),
            shell=True
    )
    subprocess.check_call(
            'gmx editconf -f "{}.gro" -o "{}.ec.gro" -n "{}" -c -d 1.0 -bt cubic'.
                format(gmx_prefix, gmx_prefix, gmx_prefix),
            shell=True
    )
